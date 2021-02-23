function plot_timeseries()
%PLOT_TIMESERIES(tg_str) plot filtered GBS time series (HDF submission)


% tg_str='NO2';

if ismac
    save_dir='/Users/raminaalwarda/Desktop/PhysicsPhD/NDACC/HDF4_data_submission/archive_v3/plots/';
elseif isunix
    save_dir='/home/kristof/work/NDACC/HDF4_data_submission/archive_v3/plots/';
end

for tg=1:3
    
    if tg==1
        tg_str='O3';
    elseif tg==2
        tg_str='NO2';
    else
        tg_str='NO2_UV';
    end
    
    if ismac
        load(['/Users/raminaalwarda/Desktop/PhysicsPhD/NDACC/HDF4_data_submission/archive_v3/PEARL-GBS_' tg_str '_VCD_all_v3.mat'])
        p=reanalysis;
        load(['/Users/raminaalwarda/Desktop/PhysicsPhD/NDACC/HDF4_data_submission/archive_v3/UT-GBS_' tg_str '_VCD_all_v3.mat'])
        ut=reanalysis;
    elseif isunix
        load(['/home/kristof/work/NDACC/HDF4_data_submission/archive_v3/PEARL-GBS_' tg_str '_VCD_all_v3.mat'])
        p=reanalysis;
        load(['/home/kristof/work/NDACC/HDF4_data_submission/archive_v3/UT-GBS_' tg_str '_VCD_all_v3.mat'])
        ut=reanalysis;
    end
        

    %% plot yearly data
    for yy=1999:2020

        ind_p=find(p.year==yy);
        ind_ut=find(ut.year==yy);

        figure(1)
        hold on

        label_count=0;
        if ~isempty(ind_ut)
            gscatter(ut.fractional_time(ind_ut),ut.mean_vcd(ind_ut),ut.ampm(ind_ut),'br','..')
            label_count=label_count+1;
        end
        if ~isempty(ind_p)
            gscatter(p.fractional_time(ind_p),p.mean_vcd(ind_p),p.ampm(ind_p),'br','oo')
            label_count=label_count+1;
        end

        if label_count==2
            legend('UT-GBS am', 'UT-GBS pm', 'PEARL-GBS am', 'PEARL-GBS pm')
        elseif label_count==1 && ~isempty(ind_ut)
            legend('UT-GBS am', 'UT-GBS pm')
        elseif label_count==1 
            legend('PEARL-GBS am', 'PEARL-GBS pm')
        else
            continue
        end

        xlim([40,310])

        switch tg_str
            case 'O3'
            ylim([0.5,1.5]*1e19)
            case 'NO2'
            ylim([0,8]*1e15)
            case 'NO2_UV'
            ylim([0,8]*1e15)
            
        end

        grid on

        xlabel(['Fractional day, ' num2str(yy) ' (UTC)'])
        ylabel([tg_str ' VCD (molec/cm^2)'])

        saveas(gcf,[save_dir tg_str '_' num2str(yy) '.png']);

        cla
        
    end

    %% plot total timeseries

    if tg~=3
        figure(2)
        subplot(2,1,tg)

        plot(p.year+p.fractional_time./daysinyear(p.year),p.mean_vcd,'bo','markersize',4), hold on
        plot(ut.year+ut.fractional_time./daysinyear(ut.year),ut.mean_vcd,'ro','markersize',4)

        ylabel([tg_str ' VCD (molec/cm^2)'])

        if tg==2
            legend('PEARL-GBS', 'UT-GBS','location','northwest')
            xlabel('Fractional year (UTC)')
            set(gcf, 'Position', [100, 100, 1200, 700]);
            saveas(gcf,[save_dir 'O3_NO2_all.png']);
        end

    else
        
        figure(2)

        plot(p.year+p.fractional_time./daysinyear(p.year),p.mean_vcd,'bo','markersize',4), hold on
        plot(ut.year+ut.fractional_time./daysinyear(ut.year),ut.mean_vcd,'ro','markersize',4)

        ylabel([tg_str ' VCD (molec/cm^2)'])

        legend('PEARL-GBS', 'UT-GBS','location','northwest')
        xlabel('Fractional year (UTC)')
        set(gcf, 'Position', [100, 100, 1200, 700]);
        saveas(gcf,[save_dir 'NO2_UV_all.png']);
        
end
end
