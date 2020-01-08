% plot GBS VCD data from a given year for reports and campaign meetings
% code assumes data to be plotted is:
%
%   UT-GBS O3 Vis on one plot
%   UT-GBS NO2 Vis, PEARL-GBS NO2 UV on another plot

%% setup

% root directory for VCD results
vcd_dir_loc='/home/kristof/work/GBS/VCD_results/';

% VCD processing state (1: final yearly data, 0: RD results)
% ordered as [UT_O3, UT_NO2, PEARL_NO2]
final_data=[0,1,1];

% year
year='2019'; % set to empty string to include all files

% save details
save_figs=1; % 0: don't save 1: save as pdf, 2: save as jpg
savedir='/home/kristof/work/campaigns/campaign_2019/report/';

%% plotting
plot_symbols={'.','o'};
plot_msize=[18,6];
no2_legend_loc='Northwest';

% loop over the two plots (o3 and no2)
for plot_num=1:2

    if plot_num==1 % ozone plot
        
        instr={'UT-GBS'};
        tgstr={'O3'};
        data_state=final_data(1);
        savename=[savedir 'ut_o3'];
        
    elseif plot_num==2 % no2 plot
        
        instr={'UT-GBS','PEARL-GBS'};
        tgstr={'NO2','NO2_UV'};
        data_state=final_data(2:3);
        savename=[savedir 'ut_p_no2'];
        
    end
    
    % set up the plot
    figure, hold on
    
    for i=1:length(instr) % for multiple datasets on the same plot
    
        % select VCD directory and file name
        if data_state(i)==1
            vcd_dir=vcd_dir_loc;
            fname=[vcd_dir instr{i} '_' tgstr{i} '_VCD_all.mat'];
        else
            vcd_dir=[vcd_dir_loc 'NDACC_RD/'];
            fname=[vcd_dir instr{i} '_' tgstr{i} '_VCD_' year 'all.mat'];
        end
        
        load(fname)
            
        reanalysis(reanalysis.year~=str2double(year),:)=[];
        if isempty(reanalysis),
            error(['Must add ' year ' ' instr{i} '_' tgstr{i} ...
                   ' data to combined dataset (use merge_reanalysis_VCDs.m)'])
        end
        
        % plot am, pm data (gscatter throws error with datetime)
        dates=mjd2k_to_date(reanalysis.mjd2k);
        ind_am=reanalysis.ampm==0;
        ind_pm=reanalysis.ampm==1;
        
        plot(dates(ind_am),reanalysis.mean_vcd(ind_am),...
             ['b' plot_symbols{i}],'markersize',plot_msize(i))
        plot(dates(ind_pm),reanalysis.mean_vcd(ind_pm),...
             ['r' plot_symbols{i}],'markersize',plot_msize(i))
        
    end
    
    % finalize plot
    grid on, box on
    
    xlabel(['Date, ' year ' (UTC)'])
    if plot_num==1
        ylabel('O_3 VCD (molec/cm^2)')
        legend('UT-GBS am','UT-GBS pm')
    else
        ylabel('NO_2 VCD (molec/cm^2)')
        legend('UT-GBS am','UT-GBS pm','PEARL-GBS am','PEARL-GBS pm',...
               'location',no2_legend_loc)
    end
    
    margin=(dates(end)-dates(1))/20;
    xlim([datenum(min(dates)-margin),datenum(max(dates)+margin)])
    
    
    set(findall(gcf,'-property','FontSize'),'FontSize',15)
    set(findall(gcf,'-property','FontName'),'FontName','Arial') 
    set(gcf, 'Position', [100, 100, 1100, 650]);
    
    % save figures
    if save_figs
        
        h=gcf;

        set(h,'Units','Inches');

        pos = get(h,'Position');

        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

        if save_figs==1
            % pdf images
            f_out=[savename '.pdf'];
            print(h,f_out,'-dpdf','-r0')
        elseif save_figs==2
            % jpg images
            f_out=[savename '.jpg'];
            print(h,f_out,'-djpeg','-r400','-opengl')
%             % png images
%             f_out=['/home/kristof/work/documents/paper_bro/figures/' fname '.png'];
%             print(h,f_out,'-dpng','-r400','-opengl')
        end
        
        close(h)
        
    end
    
    
end


