% Create profiles of (1) absolute plume melt and (2) ambient melt (where v
% = 0.01, 0.06, 0.1) from a CTD cast at the mouth of a fjord versus
% near-glacier.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here,'third_party','plume_modification')));

plume_parameters
p.plume_type = 'line';   % from Tom's workflow_example.m
p.E_0 = 0.1; % from Tom's example. Validate?
p.w   = 250; % FjordRPM auto-value. May be higher for Store as it is a large glacier.

ctd_dir = 'C:\Users\efc4\OneDrive - University of St Andrews\Desktop\PhD\First_Year\Projects\Fjord_Shelf_Observations\Working_data\OMG_data\GSW_CTDs_v3';
files  = {'OMG_Ocean_AXCTD_L2_20160915164446_gsw.nc', ...   % near glacier
          'OMG_Ocean_AXCTD_L2_20160915161939_gsw.nc'};      % fjord mouth
labels = {'Near glacier','Fjord mouth'};

d_gl   = 470;                 % grounding line depth, m (positive down)
Q      = 100;                 % total discharge, m^3/s -- Gauge from Mankoff
v_list = [0.01 0.06 0.1];     % ambient velocities, m/s

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

    zq = (0:1:floor(max(z_cast)))';                 % regular 1 m grid from the surface
    Tq = fillmissing(interp1(z_cast,T_cast,zq),'nearest');
    Sq = fillmissing(interp1(z_cast,S_cast,zq),'nearest');

    [A,~,zp] = run_plume_model(Q,Tq,Sq,-zq,0,-d_gl,p);   % plume melt, m/d in A(:,3)

    in  = zq <= d_gl;                                    % surface to grounding line
    amb = zeros(nnz(in),numel(v_list));
    for j = 1:numel(v_list)
        amb(:,j) = ambient_melt(Tq(in),Sq(in),-zq(in),v_list(j),p);
    end

    subplot(1,2,n); hold on; axis ij; box on
    plot(A(:,3),-zp,'k','LineWidth',1.5)
    plot(amb,zq(in))
    title(labels{n}); xlabel('Melt rate (m/d)'); ylabel('Depth (m)')
    legend(['Plume'; compose('Ambient v = %.2f m/s',v_list(:))],'Location','southeast')
end