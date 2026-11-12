function out = melt_profiles(fjord_id, outlet_id, opts)
% MELT_PROFILES  Plume and ambient melt profiles, near-glacier vs fjord mouth.
%   out = melt_profiles("XXX(2016)", XXXXXX)
%   out = melt_profiles("284(2016)", XXXXXX, m2019_id=XXX)    % M2019 ID used for labels only
%   out = melt_profiles("284(2016)", [80710 80711])          % several outlets, discharge summed
%   out = melt_profiles("284(2016)", 80710, ndays=30)        % name-value options (R2021a+)
%   out = melt_profiles("284(2016)", [], Q=2.5)              % manual discharge, no Mankoff lookup
% The casts table (path in config_paths) needs one row per fjord with columns
%   fjord_id, near_glacier_filename, fjord_mouth_filename, grounding_line_depth_m
% (optional: glacier_name). outlet_id is the Mankoff outlet ID (the "station" in the
% discharge file). Needs ambient_melt, load_ctd_nc, outlet_discharge_mean and
% Tom Cowton's plume_modification under third_party.
arguments
    fjord_id (1,1) string
    outlet_id double
    opts.m2019_id double = []                    % Mouginot 2019 glacier ID, label only
    opts.table (1,1) string = ""
    opts.ndays (1,1) double = 14                 % discharge window before the near-glacier cast
    opts.Q double = []                           % m^3/s; if given, replaces the Mankoff lookup
    opts.E_0 (1,1) double = 0.1                  % entrainment coefficient (Tom's example value)
    opts.v_list (1,:) double = [0.01 0.06 0.1]   % ambient velocities, m/s
    opts.w_list (1,:) double = [50 100 250]      % plume widths, m (plume depends only on Q/w)
    opts.save_figs (1,1) logical = true
end
P = config_paths();
if opts.table == "", opts.table = string(P.table); end

%% ---- look up this fjord in the table ----
T   = readtable(opts.table,'TextType','string');
row = T(T.fjord_id == fjord_id,:);
assert(height(row)==1, 'fjord_id "%s" matched %d rows in %s', fjord_id, height(row), opts.table)
files = {char(row.near_glacier_filename), char(row.fjord_mouth_filename)};
d_gl  = row.grounding_line_depth_m;
assert(~isnan(d_gl), 'No grounding line depth for %s', fjord_id)

name = char(fjord_id);
if ismember('glacier_name',T.Properties.VariableNames), name = sprintf('%s (%s)',row.glacier_name,fjord_id); end
if ~isempty(opts.m2019_id), name = sprintf('%s, M2019 ID %d',name,opts.m2019_id); end

%% ---- discharge ----
if ~isempty(opts.Q)
    Q = opts.Q;
    qtag = 'Qmanual';
else
    assert(~isempty(outlet_id), 'Give the Mankoff outlet_id, or pass Q=...')
    tok  = regexp(files{1},'_(\d{14})_','tokens','once');          % cast time from filename
    f_nc = fullfile(P.mankoff_dir, sprintf('MAR_%s.nc',tok{1}(1:4)));
    Q = outlet_discharge_mean(f_nc, outlet_id, files{1}, opts.ndays);
    qtag = sprintf('%dd',opts.ndays);
end
fprintf('%s: d_gl = %d m, Q = %.2f m^3/s\n', name, d_gl, Q)

%% ---- model parameters ----
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here,'third_party','plume_modification')));
p = struct();
plume_parameters
p.plume_type = 'line';
p.E_0 = opts.E_0;
p.w   = opts.w_list(1);          % placeholder; overridden per run below

%% ---- setup ----
labels = {'Near glacier','Fjord mouth'};
ylims  = [0 d_gl];
v_list = opts.v_list;  w_list = opts.w_list;
lsty   = {'-','--',':','-.'};
zc     = (0:1:d_gl)';
res    = struct('plume',{[],[]},'amb',{[],[]},'TF',{[],[]});
tag    = sprintf('%s_%s',regexprep(char(fjord_id),'[^\w]','_'),qtag);

%% ---- figure 1: melt profiles for each cast ----
fig1 = figure;
for n = 1:2
    [z_cast,T_cast,S_cast] = load_ctd_nc(fullfile(P.ctd_dir,files{n}), ...
        'depth','conservative_temperature','absolute_salinity');
    
    fprintf('%s: depth %.1f to %.1f m, %d points\n', ...
        labels{n},min(z_cast),max(z_cast),numel(z_cast))

    % Check whether CTD reaches the grounding line
    if max(z_cast) < d_gl
        warning(['%s cast reaches %.1f m, short of the %.1f m ' ...
            'grounding line. Extrapolating deepest T/S values.'], ...
            labels{n},max(z_cast),d_gl)
    end

    % Interpolate CTD profiles onto a 1 m depth grid
    zq = (0:1:ceil(d_gl))';

    Tq = interp1(z_cast,T_cast,zq,'linear','extrap');
    Sq = interp1(z_cast,S_cast,zq,'linear','extrap');

    % Use constant bottom values rather than linear extrapolation
    % below the deepest CTD measurement
    Tq(zq > max(z_cast)) = T_cast(end);
    Sq(zq > max(z_cast)) = S_cast(end);

    in  = zq <= d_gl;
    amb = zeros(nnz(in),numel(v_list));
    for j = 1:numel(v_list)
        amb(:,j) = ambient_melt(Tq(in),Sq(in),-zq(in),v_list(j),p);
    end

    subplot(1,2,n); hold on; axis ij; box on
    h_amb = plot(amb,zq(in));

    h_pl = gobjects(1,numel(w_list));
    res(n).plume = nan(numel(zc),numel(w_list));
    for k = 1:numel(w_list)
        pw = p;  pw.w = w_list(k);
        [A,~,zp] = run_plume_model(Q,Tq,Sq,-zq,0,-d_gl,pw);   % A(:,3) = melt, m/d
        h_pl(k) = plot(A(:,3),-zp,'k','LineStyle',lsty{mod(k-1,numel(lsty))+1},'LineWidth',1.5);
        [zpd,iu] = unique(-zp);  mp = A(:,3);
        res(n).plume(:,k) = interp1(zpd,mp(iu),zc);
        fprintf('  w = %3d m: plume ends at %.0f m, max melt %.2f m/d\n',w_list(k),-zp(end),max(A(:,3)))
    end
    res(n).amb = interp1(zq(in),amb,zc);
    res(n).TF  = interp1(zq,Tq-(p.lambda1*Sq+p.lambda2+p.lambda3*(-zq)),zc);

    title(labels{n}); xlabel('Melt rate (m/d)'); ylabel('Depth (m)')
    legend([h_pl h_amb'],[compose('Plume, w = %d m',w_list(:))' ...
                          compose('Ambient, v = %.2f m/s',v_list(:))'],'Location','southeast')
end
linkaxes(findobj(fig1,'Type','axes'),'xy'); ylim(ylims)
sgtitle(sprintf('%s (Q = %.2f m^3/s)',name,Q))
if opts.save_figs, exportgraphics(fig1,sprintf('melt_profiles_%s.png',tag),'Resolution',300); end

%% ---- figure 2: differences, fjord mouth minus near glacier ----
if ~isempty(res(1).amb) && ~isempty(res(2).amb)
    fig2 = figure;
    subplot(1,3,1); hold on; axis ij; box on
    dpl = res(2).plume - res(1).plume;
    for k = 1:numel(w_list)
        plot(dpl(:,k),zc,'k','LineStyle',lsty{mod(k-1,numel(lsty))+1},'LineWidth',1.5)
    end
    xline(0,'k:'); xlabel('\Delta plume melt (m/d)'); ylabel('Depth (m)'); title('Plume')
    legend(compose('w = %d m',w_list(:)),'Location','southeast')

    subplot(1,3,2); hold on; axis ij; box on
    plot(res(2).amb-res(1).amb,zc); xline(0,'k:')
    xlabel('\Delta ambient melt (m/d)'); title('Ambient')
    legend(compose('v = %.2f m/s',v_list(:)),'Location','southeast')

    subplot(1,3,3); hold on; axis ij; box on
    plot(res(2).TF-res(1).TF,zc,'r'); xline(0,'k:')
    xlabel('\Delta thermal forcing (\circC)'); title('Thermal forcing')

    sgtitle([name ': fjord mouth minus near glacier (positive = higher at mouth)'])
    linkaxes(findobj(fig2,'Type','axes'),'y'); ylim(ylims)
    if opts.save_figs, exportgraphics(fig2,sprintf('melt_differences_%s.png',tag),'Resolution',300); end
end

out = struct('fjord_id',fjord_id,'outlet_id',outlet_id,'m2019_id',opts.m2019_id,'name',name, ...
             'Q',Q,'d_gl',d_gl,'zc',zc,'res',res,'v_list',v_list,'w_list',w_list,'opts',opts);
end

% CALL
% out = melt_profiles("<fjord_id>", "<Mankoff_gate_id>", m2019_id="<insert>");
