% check_ctds.m  Quick look at the two CTD casts before running the melt model
ctd_dir = 'C:\Users\efc4\OneDrive - University of St Andrews\Desktop\PhD\First_Year\Projects\Fjord_Shelf_Observations\Working_data\OMG_data\GSW_CTDs_v3';
files  = {'OMG_Ocean_AXCTD_L2_20160915164446_gsw.nc', ...
    'OMG_Ocean_AXCTD_L2_20160915161939_gsw.nc'};
labels = {'Near glacier','Fjord mouth'};
cols   = {'r','b'};
d_gl   = 470;

figure
for n = 1:2
    [z,T,S] = load_ctd_nc(fullfile(ctd_dir,files{n}), ...
        'depth','conservative_temperature','absolute_salinity');

    zq = (0:1:floor(max(z)))';
    Tq = fillmissing(interp1(z,T,zq),'nearest');
    Sq = fillmissing(interp1(z,S,zq),'nearest');

    fprintf('%s: T %.2f to %.2f, S %.2f to %.2f, first depth %.2f, NaNs before fill: T %d, S %d\n', ...
        labels{n},min(T),max(T),min(S),max(S),min(z),sum(isnan(interp1(z,T,zq))),sum(isnan(interp1(z,S,zq))))

    subplot(1,3,1); hold on; axis ij; box on
    plot(T,z,'.','Color',cols{n},'MarkerSize',2,'HandleVisibility','off')
    plot(Tq,zq,'-','Color',cols{n},'DisplayName',labels{n})
    xlabel('Conservative temperature (\circC)'); ylabel('Depth (m)')

    subplot(1,3,2); hold on; axis ij; box on
    plot(S,z,'.','Color',cols{n},'MarkerSize',2,'HandleVisibility','off')
    plot(Sq,zq,'-','Color',cols{n},'DisplayName',labels{n})
    xlabel('Absolute salinity (g/kg)')

    subplot(1,3,3); hold on; axis ij; box on
    plot(gsw_sigma0(Sq,Tq),zq,'-','Color',cols{n},'DisplayName',labels{n})
    xlabel('\sigma_0 (kg/m^3)')
end
for k = 1:3
    subplot(1,3,k); yline(d_gl,'k--','Grounding line','HandleVisibility','off');
end
subplot(1,3,1); legend('Location','southeast')