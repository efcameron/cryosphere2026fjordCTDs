% Create profiles of (1) absolute plume melt and (2) ambient melt (where v
% = 0.01, 0.06, 0.1) from a CTD cast at the mouth of a fjord versus
% near-glacier.

%% ---- paths and parameters ----
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here,'third_party','plume_modification')));

plume_parameters
p.plume_type = 'line';   % from Tom's workflow_example.m
p.E_0 = 0.1;             % from Tom's example. Validate?
p.w   = 250;             % placeholder only: overridden by w_list below, but the field must exist

%% ---- inputs ----
ctd_dir = 'C:\Users\efc4\OneDrive - University of St Andrews\Desktop\PhD\First_Year\Projects\Fjord_Shelf_Observations\Working_data\OMG_data\GSW_CTDs_v3';
files  = {'OMG_Ocean_AXCTD_L2_20160915164446_gsw.nc', ...   % near glacier
          'OMG_Ocean_AXCTD_L2_20160915161939_gsw.nc'};      % fjord mouth
labels = {'Near glacier','Fjord mouth'};

d_gl   = 470;                 % grounding line depth, m (positive down)
ylims  = [0 d_gl];  
Q      = 2.5;                 % m^3/s, 30-day mean before the CTDs (Mankoff 2020); 14-day mean is ~2
v_list = [0.01 0.06 0.1];     % ambient velocities, m/s
w_list = [50 100 250];        % plume widths, m (the plume depends only on Q/w)
lsty   = {'-','--',':'};      % line styles for the three plume widths

zc  = (0:1:d_gl)';            % common depth grid for differencing
res = struct('plume',{[],[]},'amb',{[],[]},'TF',{[],[]});

%% ---- figure 1: melt profiles for each cast ----
figure
for n = 1:2
    [z_cast,T_cast,S_cast] = load_ctd_nc(fullfile(ctd_dir,files{n}), ...
        'depth','conservative_temperature','absolute_salinity');
    fprintf('%s: depth %.1f to %.1f m, %d points\n',labels{n}, ...
        min(z_cast),max(z_cast),numel(z_cast))
    if max(z_cast) < d_gl
        warning('%s cast stops above the grounding line (%.0f m).',labels{n},max(z_cast))
        continue
    end

    % regrid the cast to a regular 1 m grid from the surface
    zq = (0:1:floor(max(z_cast)))';
    Tq = fillmissing(interp1(z_cast,T_cast,zq),'nearest');
    Sq = fillmissing(interp1(z_cast,S_cast,zq),'nearest');

    % ambient melt, one curve per velocity
    in  = zq <= d_gl;
    amb = zeros(nnz(in),numel(v_list));
    for j = 1:numel(v_list)
        amb(:,j) = ambient_melt(Tq(in),Sq(in),-zq(in),v_list(j),p);
    end

    subplot(1,2,n); hold on; axis ij; box on
    h_amb = plot(amb,zq(in));

    % plume melt, one curve per width
    h_pl = gobjects(1,numel(w_list));
    res(n).plume = nan(numel(zc),numel(w_list));
    for k = 1:numel(w_list)
        pw = p;  pw.w = w_list(k);                           % copy p, change only the width
        [A,~,zp] = run_plume_model(Q,Tq,Sq,-zq,0,-d_gl,pw); % A(:,3) = melt, m/d
        h_pl(k) = plot(A(:,3),-zp,'k','LineStyle',lsty{k},'LineWidth',1.5);
        [zpd,iu] = unique(-zp);  mp = A(:,3);
        res(n).plume(:,k) = interp1(zpd,mp(iu),zc);          % NaN where the plume never reaches
        fprintf('  w = %3d m: plume ends at %.0f m, max melt %.2f m/d\n', ...
            w_list(k),-zp(end),max(A(:,3)))
    end

    % store on the common grid for the difference plots
    res(n).amb = interp1(zq(in),amb,zc);                     % columns = v_list
    res(n).TF  = interp1(zq,Tq-(p.lambda1*Sq+p.lambda2+p.lambda3*(-zq)),zc);

    title(labels{n}); xlabel('Melt rate (m/d)'); ylabel('Depth (m)')
    legend([h_pl h_amb'], [compose('Plume, w = %d m',w_list(:))' ...
                           compose('Ambient, v = %.2f m/s',v_list(:))'], ...
           'Location','southeast')
end
linkaxes(findobj(gcf,'Type','axes'),'xy') % same scales on both casts
ylim(ylims)

%% ---- figure 2: differences, fjord mouth minus near glacier ----
if ~isempty(res(1).amb) && ~isempty(res(2).amb)
    figure

    subplot(1,3,1); hold on; axis ij; box on
    dpl = res(2).plume - res(1).plume;
    for k = 1:numel(w_list)
        plot(dpl(:,k),zc,'k','LineStyle',lsty{k},'LineWidth',1.5)
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

    sgtitle('Fjord mouth minus near glacier (positive = higher at mouth)')
    linkaxes(findobj(gcf,'Type','axes'),'y')                 % share depth axis across the three panels
    ylim(ylims)
end