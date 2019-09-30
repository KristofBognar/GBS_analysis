% plot reanalysis VCDs and compare to Cristen's data

tg=2;

savedir='/home/kristof/work/GBS/VCD_results/plots/';
load('/home/kristof/work/GBS/VCD_results/VCDs_Cristen_JUNE2011.mat');

if tg==1
    du=2.687e16;

    load('/home/kristof/work/GBS/VCD_results/UT-GBS_O3_VCD_all.mat');
    utgbs=reanalysis;

    load('/home/kristof/work/GBS/VCD_results/VCD_with_bad_rcd/PEARL-GBS_O3_VCD_all.mat');
    pgbs=reanalysis;

    %% plot all
    
    % UT-GBS
    time=utgbs.year+((utgbs.fd-1)./365);
    time_c=vcd_u1_o3_ndacc_filt(:,1)+((vcd_u1_o3_ndacc_filt(:,4)-1)./365);
    
    gscatter(time,utgbs.mean_vcd./du,utgbs.ampm,'br','..',18); hold on
    gscatter(time_c,vcd_u1_o3_ndacc_filt(:,13)./du, vcd_u1_o3_ndacc_filt(:,3), 'br', 'oo')
    
    xlabel('Fractional year (UTC)')
    ylabel('Ozone VCD (DU)')
    xlim([1999,2017])
    legend('UT-GBS am','UT-GBS pm','Cristen am','Cristen pm')

    savename= ['O3_VCD_UT-GBS_vs_Cristen.png'];
    f=gcf;
    print(f,fullfile(savedir,savename),'-dpng','-r300','-opengl') %save file
    savefig([savedir 'O3_VCD_UT-GBS_vs_Cristen.fig'])
    clf
 
    % PEARL-GBS
    time=pgbs.year+((pgbs.fd-1)./365);
    time_c=vcd_p1_o3_ndacc_filt(:,1)+((vcd_p1_o3_ndacc_filt(:,4)-1)./365);
    
    gscatter(time,pgbs.mean_vcd./du,pgbs.ampm,'br','..',18); hold on
    gscatter(time_c,vcd_p1_o3_ndacc_filt(:,13)./du, vcd_p1_o3_ndacc_filt(:,3), 'br', 'oo')
    
    xlabel('Fractional year (UTC)')
    ylabel('Ozone VCD (DU)')
    xlim([2006,2010])
    legend('PEARL-GBS am','PEARL-GBS pm','Cristen am','Cristen pm')

    savename= ['O3_VCD_PEARL-GBS_vs_Cristen.png'];
    f=gcf;
    print(f,fullfile(savedir,savename),'-dpng','-r300','-opengl') %save file
    savefig([savedir 'O3_VCD_PEARL-GBS_vs_Cristen.fig'])
    clf
    
    %% plot yearly

    for i=1999:2016

        if i==2001 || i==2002, continue, end
        
        ind_u=find(utgbs.year==i);
        ind_p=find(pgbs.year==i);

        gscatter(utgbs.fd(ind_u),utgbs.mean_vcd(ind_u)./du,utgbs.ampm(ind_u),'br','..',18); hold on
        if ~isempty(ind_p)
            gscatter(pgbs.fd(ind_p),pgbs.mean_vcd(ind_p)./du,pgbs.ampm(ind_p),'br','oo');
        end

        xlabel(['Fractional day, ' num2str(i) ' (UTC)'])
        ylabel('Ozone VCD (DU)')
        if isempty(ind_p)
            legend('UT-GBS am','UT-GBS pm','location','best')
        else
            legend('UT-GBS am','UT-GBS pm','PEARL-GBS am','PEARL-GBS pm','location','best')
        end

        savename= ['O3_VCD_' num2str(i) '.png'];
        f=gcf;
        print(f,fullfile(savedir,savename),'-dpng','-r300','-opengl') %save file
        clf
        
    end

elseif tg==2
    du=1;

    load('/home/kristof/work/GBS/VCD_results/UT-GBS_NO2_VCD_all.mat');
    utgbs=reanalysis;

    load('/home/kristof/work/GBS/VCD_results/VCD_with_bad_rcd/PEARL-GBS_NO2_VCD_all.mat');
    pgbs=reanalysis;

    %% plot all

    % UT-GBS
    time=utgbs.year+((utgbs.fd-1)./365);
    time_c=vcd_u1_no2_ndacc_filt(:,1)+((vcd_u1_no2_ndacc_filt(:,4)-1)./365);
    
    gscatter(time,utgbs.mean_vcd./du,utgbs.ampm,'br','..',18); hold on
    gscatter(time_c,vcd_u1_no2_ndacc_filt(:,13)./du, vcd_u1_no2_ndacc_filt(:,3), 'br', 'oo')
    
    xlabel('Fractional year (UTC)')
    ylabel('NO_2 VCD (molec/cm^2)')
    xlim([1999,2017])
    legend('UT-GBS am','UT-GBS pm','Cristen am','Cristen pm')

    savename= ['NO2_VCD_UT-GBS_vs_Cristen.png'];
    f=gcf;
    print(f,fullfile(savedir,savename),'-dpng','-r300','-opengl') %save file
    savefig([savedir 'NO2_VCD_UT-GBS_vs_Cristen.fig'])
    clf
 
    % PEARL-GBS
    time=pgbs.year+((pgbs.fd-1)./365);
    time_c=vcd_p1_no2_ndacc_filt(:,1)+((vcd_p1_no2_ndacc_filt(:,4)-1)./365);
    
    gscatter(time,pgbs.mean_vcd./du,pgbs.ampm,'br','..',18); hold on
    gscatter(time_c,vcd_p1_no2_ndacc_filt(:,13)./du, vcd_p1_no2_ndacc_filt(:,3), 'br', 'oo')
    
    xlabel('Fractional year (UTC)')
    ylabel('NO_2 VCD (molec/cm^2)')
    xlim([2006,2010])
    legend('PEARL-GBS am','PEARL-GBS pm','Cristen am','Cristen pm')

    savename= ['NO2_VCD_PEARL-GBS_vs_Cristen.png'];
    f=gcf;
    print(f,fullfile(savedir,savename),'-dpng','-r300','-opengl') %save file
    savefig([savedir 'NO2_VCD_PEARL-GBS_vs_Cristen.fig'])
    clf

    %% plot yearly

    for i=1999:2016

        if i==2001 || i==2002, continue, end
        
        ind_u=find(utgbs.year==i);
        ind_p=find(pgbs.year==i);

        gscatter(utgbs.fd(ind_u),utgbs.mean_vcd(ind_u)./du,utgbs.ampm(ind_u),'br','..',18); hold on
        if ~isempty(ind_p)
            gscatter(pgbs.fd(ind_p),pgbs.mean_vcd(ind_p)./du,pgbs.ampm(ind_p),'br','oo');
        end

        xlabel(['Fractional day, ' num2str(i) ' (UTC)'])
        ylabel('NO_2 VCD (molec/cm^2)')
        if isempty(ind_p)
            legend('UT-GBS am','UT-GBS pm','location','best')
        else
            legend('UT-GBS am','UT-GBS pm','PEARL-GBS am','PEARL-GBS pm','location','best')
        end

        savename= ['NO2_VCD_' num2str(i) '.png'];
        f=gcf;
        print(f,fullfile(savedir,savename),'-dpng','-r300','-opengl') %save file
        clf
        
    end

end

clearvars

