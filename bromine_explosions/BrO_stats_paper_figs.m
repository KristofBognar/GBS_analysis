function BrO_stats_paper_figs()
%BRO_STATS_PAPER_FIGS to plot figures for paper describing BrO statistics
%in Eureka
%
% This code contains modified versions of the code from several other
% functions -- collected here to easily update/reproduce figures in the
% paper.
%
% Some data used in plotting are not the usual output files, but copies
% saved in /home/kristof/work/documents/paper_bro 


%% setup

% what to plot
presentation_plots=0;

windrose=0;
plot_box=0;
plot_box_weather=0;
bro_dailymean=0;
bro_dailymean_o3=1;
plot_pca=0;
weather_corr=0;
o3_aer_wspd=0;
sens_map=0;
plot_ssa=0;

si_contact_log=0;
plot_log_si=1; % log scale for SI contact axes
btraj_len='3'; % length of back trajectories

% 0 to not save, 1 to save as pdf, 2 to save as jpg
save_figs=2;

% uniform look
fig_fs=14; % font size on figures
arr_txt_fs=9;
fig_font='Arial'; % font for plotting

% line with for box plot, and markers
box_lw=1.5;
box_outlier=4; % outlier size for box plots
    
% plotting limits for BrO scatter plots 
bro_lim_full=16.5; %*e13, all data
bro_lim=12.5; %*e13, 7 values are above 8e13
%%% v1 plots with all 1ppt a priori surf conc: bro_lim_full=12.5; bro_lim=8;


% plotting aids
% c_list={'r','g','b','m','c'};
% % c_list2=[[1 0 0];...
% %          [0 1 0];...
% %          [0 0 1];...
% %          [0 1 1]];
% c_list2=[[0 0 0.4];...
%          [0 0 1];...
%          [0 0.6 1];...
%          [0 1 0.5]];

c_scale='winter';     
c_list=parula(5);     
c_list=c_list(1:4,:);     

m_names={'March','April','May'};
letters='abcdefghijklmnopqrstuvwxyz';


%% load and filter BEE dataset
load('/home/kristof/work/BEEs/BEE_dataset_all.mat')

% exclude 2015 (intermittent measurements)
bee_dataset(bee_dataset.times.Year==2015,:)=[];

%%% test filters, not used in paper
% optimistic detction limit
% bee_dataset(bee_dataset.bro_col-bee_dataset.bro_col_err*3<=0,:)=[];
% pessimistic detection limit
% bee_dataset(bee_dataset.bro_col<mean(bee_dataset.bro_col_err)*3,:)=[];
% 'high' BrO only 
% bee_dataset(bee_dataset.bro_col<prctile(bee_dataset.bro_col,50),:)=[]; % above Q3
% bee_dataset(bee_dataset.wspd_ms<8,:)=[];
% bee_dataset(bee_dataset.o3_surf<=15,:)=[];
% bee_dataset(bee_dataset.length_3day<=750,:)=[];
% bee_dataset(bee_dataset.mixing_height_3day>300,:)=[];
%%%

% setup wdir indices
% could fill NaNs with EWS data, but actual wdir is different from PWS ~70%
% of the time (even though the wind rose looks similar)
ind_N=bee_dataset.N_SE_rest==1;
ind_SE=bee_dataset.N_SE_rest==2;
ind_rest=bee_dataset.N_SE_rest==3;

% load data grouped for flexpart
load('/home/kristof/work/BEEs/BEE_dataset_flexpart.mat')

% filter dataset later (sensitivities are saved in a way that all runs must
% be included in the indexing)



% apriori=get_BrO_apriori(bee_dataset.times);
% bee_dataset.bro_col(apriori.surf_ppt==5)=bee_dataset.bro_col(apriori.surf_ppt==5)*0.8;

%% windrose data
if windrose

    % PWS wind data
    load('/home/kristof/work/weather_stations/ridge_lab/PWS_all.mat');
    data(data.DateTime.Year==2015 | (month(data.DateTime)>5 | month(data.DateTime)<3),:)=[];

    WindRose(data.WindDir,data.WindSpd,'AngleNorth',0,'AngleEast',90,...
             'FreqLabelAngle',45,'vWinds',[0 3 6 9 12 15],'MaxFrequency',15,...
             'height',250,'width',700,'LabLegend','');
    
    set(gcf, 'Position', [100, 100, 750, 500]);
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)

    % windrose polygons are messed up in PDF
    if save_figs, save_pdf(2, 'windrose'), end
         
end

%% 2x2 box plot with BrO, AOD, aer and O3
if plot_box
    
    labelx=0.78;
    labely=0.93;

    % set up figure
    figure
    set(gcf, 'Position', [100, 100, 1000, 600]);
    fig_ax = tight_subplot(2,2,[0.07,0.07],[0.08,0.1],[0.08,0.03]);
     
    plot_data=bee_dataset;
    
    bad_ind=find(bee_dataset.times.Year==2017 & bee_dataset.times.Month==4);
    plot_data.bro_col(bad_ind)=NaN;
    plot_data.aer_ext(bad_ind)=NaN;
    plot_data.aer_halfmicron(bad_ind)=NaN;
    plot_data.SMPS_100_500(bad_ind)=NaN;
    plot_data.o3_surf(bad_ind)=NaN;
    
    plot_data=[plot_data;plot_data(end,:)];
    plot_data.times(end)=datetime(2017,05,13);
    plot_data.bro_col(end)=NaN;
    plot_data.aer_ext(end)=NaN;
    plot_data.aer_halfmicron(end)=NaN;
    plot_data.SMPS_100_500(end)=NaN;
    plot_data.o3_surf(end)=NaN;
    
    box_group=plot_data.times.Year+plot_data.times.Month*10000;

    %%%%%%%%%% BrO box plot
    axes(fig_ax(1))
    
    % set positions to separate bar groups by wdir
    a=0.2; % controls width of each group
    nbars=4; % number of bars in each group 
    
    ll=plot_yearly_box(plot_data.bro_col,box_group,a,nbars,c_list,box_lw,box_outlier,1);
    
    ll.Position(1)=0.36;
    ll.Position(2)=0.92;
               
    ylabel('BrO part. col. (molec cm^{-2})')
    ylim([-0.5,bro_lim_full]*1e13)
    
    set(gca,'YtickLabel',[0,5,10,15])
    
    text(0.847,1.75e14,'\times10^{13}', 'color','k');
    text(labelx,labely,'a)','color','k','FontWeight','bold','Units','normalized')

    
    %%%%%%%%%% AOD box plot
    axes(fig_ax(2))
    
    plot_yearly_box(plot_data.aer_ext,box_group,a,nbars,c_list,box_lw,box_outlier,0);
    
    set(gca, 'YScale', 'log')
    ylb=ylabel('AOD');
    ylb.Position(1)=ylb.Position(1)+0.17;
    
    text(labelx,labely,'b)','color','k','FontWeight','bold','Units','normalized')
    
    %%%%%%%%%% ozone box plot
    axes(fig_ax(3))
    
    plot_yearly_box(plot_data.o3_surf,box_group,a,nbars,c_list,box_lw,box_outlier,0);
    
    ylabel('O_3 (ppbv)')
    
    text(labelx,labely,'c)','color','k','FontWeight','bold','Units','normalized')
    
    %%%%%%%%%% SSA box plot
    axes(fig_ax(4))
    
    plot_yearly_box(plot_data.aer_halfmicron,box_group,a,nbars,c_list,box_lw,box_outlier,0);
    ylabel('d_p > 0.5 \mum (cm^{-3})')
    
%     plot_yearly_box(plot_data.SMPS_100_500,box_group,a,nbars,c_list,box_lw,box_outlier,0);
%     ylabel('0.1 < d_p < 0.5 \mum (cm^{-3})')
    
    text(labelx,labely,'d)','color','k','FontWeight','bold','Units','normalized')    
    
    %%% plot formatting
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
%     ll.FontSize=11;

    save_pdf(save_figs, 'data_stats')
    
    
    
end

%% 2x2 box plot with weather stuff
if plot_box_weather
    
    labelx=0.78;
    labely=0.93;

    % set up figure
    figure
    set(gcf, 'Position', [100, 100, 1000, 600]);
    fig_ax = tight_subplot(2,2,[0.07,0.07],[0.08,0.1],[0.08,0.03]);
        
    plot_data=bee_dataset;
    
    bad_ind=find(bee_dataset.times.Year==2017 & bee_dataset.times.Month==4);
    plot_data.T_PWS(bad_ind)=NaN;
    plot_data.T_EWS(bad_ind)=NaN;
    plot_data.sonde_dT(bad_ind)=NaN;
    plot_data.wspd_ms(bad_ind)=NaN;
    
    plot_data=[plot_data;plot_data(end,:)];
    plot_data.times(end)=datetime(2017,05,13);
    plot_data.T_PWS(end)=NaN;
    plot_data.T_EWS(end)=NaN;
    plot_data.sonde_dT(end)=NaN;
    plot_data.wspd_ms(end)=NaN;
    
    box_group=plot_data.times.Year+plot_data.times.Month*10000;

    % set positions to separate bar groups 
    a=0.2; % controls width of each group
    nbars=4; % number of bars in each group 
    
    %%%%%%%%%% RL T box plot
    axes(fig_ax(1))
    
    ll=plot_yearly_box(plot_data.T_PWS,box_group,a,nbars,c_list,box_lw,box_outlier,1);
    
    ll.Position(1)=0.118;
    ll.Position(2)=0.92;
               
    ylabel('Temperature (°C)')

    text(labelx,labely,'a)','color','k','FontWeight','bold','Units','normalized')

    axes(fig_ax(2))
    delete(gca)
    
    %%%%%%%%%% dT
    axes(fig_ax(3))
    
    plot_yearly_box(plot_data.sonde_dT,box_group,a,nbars,c_list,box_lw,box_outlier,0);
    
    ylabel('Inversion strength (°C)')
    
    text(labelx,labely,'c)','color','k','FontWeight','bold','Units','normalized')
    
    %%%%%%%%%% SSA box plot
    axes(fig_ax(4))
    
    plot_yearly_box(plot_data.wspd_ms,box_group,a,nbars,c_list,box_lw,box_outlier,0);
    
    ylabel('Wind speed (m s^{-1})')
    
    text(labelx,labely,'d)','color','k','FontWeight','bold','Units','normalized')    
    
    %%% plot formatting
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
%     ll.FontSize=11;

    save_pdf(save_figs, 'weather_stats')
    
end

%% daily mean BrO, with min-max range
if bro_dailymean
    
%     % check if daily mean BrO was saved already
%     if ~exist('/home/kristof/work/documents/paper_bro/data/daily_mean_BrO.mat','file')

        % find unique days
        unique_days=unique(...
            [bee_dataset.times.Year,bee_dataset.times.Month,bee_dataset.times.Day],'rows');

        dmean_bro=NaN(size(unique_days,1),4);

        % loop over unique days
        for i=1:size(unique_days,1)

            tmp=bee_dataset.bro_col( bee_dataset.times.Year==unique_days(i,1) & ...
                                     bee_dataset.times.Month==unique_days(i,2) & ...
                                     bee_dataset.times.Day==unique_days(i,3)          );

            dmean_bro(i,1)=mean(tmp);
            dmean_bro(i,2)=std(tmp);
            dmean_bro(i,3)=min(tmp);
            dmean_bro(i,4)=max(tmp);

        end

        % save results in a table
        dmean_bro=array2table(dmean_bro,'variablenames',{'mean','std','min','max'}); 
        dmean_bro.DateTime=datetime(unique_days)+hours(12);
        
%         save('/home/kristof/work/documents/paper_bro/data/daily_mean_BrO.mat','dmean_bro')
%     else
%         % if already saved 
%         load('/home/kristof/work/documents/paper_bro/data/daily_mean_BrO.mat')
%     end
    
    % plot figure
    figure
%     set(gcf, 'Position', [100, 100, 1000, 600]);
%     fig_ax = tight_subplot(4,1,[0.065,0.04],[0.12,0.068],[0.088,0.04]);
    set(gcf, 'Position', [100, 100, 1000, 700]);
    fig_ax = tight_subplot(4,1,[0.065,0.04],[0.104,0.06],[0.088,0.04]);
    
    box_w=0.7; % box width in days
    
    for yr=2016:2019
    
        axes(fig_ax(yr-2015))
        
        ind=(dmean_bro.DateTime.Year==yr);

        plot_x=datenum(dmean_bro.DateTime(ind));
        plot_mean=dmean_bro.mean(ind);
        plot_min=dmean_bro.min(ind);
        plot_max=dmean_bro.max(ind);
        plot(1,1), hold on % need to have something plotted before turning hold on,
                           % tight_subplot removes ylabels otherwise...

        for dy=1:length(plot_x)
            
            % draw rectangle for each day, with min and max BrO columns as
            % the top and bottom edges
            % position is [x,y,w,h] (bottom left corner)
            rectangle('position',...
                      [plot_x(dy)-box_w/2,plot_min(dy),box_w,plot_max(dy)-plot_min(dy)],...
                      'facecolor',c_list(2,:),'edgecolor',c_list(2,:))
            rectangle('position',...
                      [plot_x(dy)-box_w/2,plot_mean(dy),box_w,1e12],...
                      'facecolor','k','edgecolor','k')
                  
        end

        grid on;
        ylim([-0.5,bro_lim_full]*1e13)
        xlim([datenum(yr,3,1),datenum(yr,6,2)])
        
%         plot(datenum([[yr,3,1];[yr,6,1]]),[3.7,3.7]*1e12)
        
        set(gca, 'XTick', datenum([[yr,3,1];[yr,3,15];[yr,4,1];[yr,4,15];...
                                   [yr,5,1];[yr,5,15];[yr,6,1]]))
        set(gca, 'XTicklabel', {'01/03','15/03','01/04','15/04','01/05','15/05','01/06'})
        
        if yr~=2019
            set(gca,'xticklabel',[])
        else
            xlabel('Date (dd/mm, UTC)')
        end
        
        if yr==2017
            ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
            ylb.Position(1)=ylb.Position(1)-1;
            ylb.Position(2)=ylb.Position(2)-max(dmean_bro.mean);
        end
        
        text(0.9,0.8,num2str(yr), 'color','k','Units','normalized',...
        'fontsize',fig_fs,'fontweight','bold')
        
    end
    
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)

    save_pdf(save_figs, 'BrO_dmean')
    
end

if bro_dailymean_o3
    
        % find unique days
        unique_days=unique(...
            [bee_dataset.times.Year,bee_dataset.times.Month,bee_dataset.times.Day],'rows');

        dmean_bro=NaN(size(unique_days,1),4);
        dmean_o3=NaN(size(unique_days,1),4);

        % loop over unique days
        for i=1:size(unique_days,1)

            tmp=bee_dataset.bro_col( bee_dataset.times.Year==unique_days(i,1) & ...
                                     bee_dataset.times.Month==unique_days(i,2) & ...
                                     bee_dataset.times.Day==unique_days(i,3)          );

            dmean_bro(i,1)=mean(tmp);
            dmean_bro(i,2)=std(tmp);
            dmean_bro(i,3)=min(tmp);
            dmean_bro(i,4)=max(tmp);

            tmp=bee_dataset.o3_surf( bee_dataset.times.Year==unique_days(i,1) & ...
                                     bee_dataset.times.Month==unique_days(i,2) & ...
                                     bee_dataset.times.Day==unique_days(i,3)          );
            dmean_o3(i,1)=mean(tmp);
            dmean_o3(i,2)=std(tmp);
            dmean_o3(i,3)=min(tmp);
            dmean_o3(i,4)=max(tmp);
            
        end

        % save results in a table
        dmean_bro=array2table(dmean_bro,'variablenames',{'mean','std','min','max'}); 
        dmean_bro.DateTime=datetime(unique_days)+hours(12);
        
        dmean_o3=dmean_o3.*0.5e13; % scale to BrO
        dmean_o3=array2table(dmean_o3,'variablenames',{'mean','std','min','max'}); 
        dmean_o3.DateTime=datetime(unique_days)+hours(12);
        
    
    % plot figure
    figure
    set(gcf, 'Position', [100, 100, 1000, 700]);
    fig_ax = tight_subplot(4,1,[0.065,0.04],[0.104,0.06],[0.088,0.1]);
    
    box_w=0.7; % box width in days
    
    for yr=2016:2019
    
        axes(fig_ax(yr-2015))
        
        ind=(dmean_bro.DateTime.Year==yr);

        plot_x=datenum(dmean_bro.DateTime(ind));
        plot_mean=dmean_bro.mean(ind);
        plot_min=dmean_bro.min(ind);
        plot_max=dmean_bro.max(ind);

        plot_mean_o3=dmean_o3.mean(ind);
        plot_min_o3=dmean_o3.min(ind);
        plot_max_o3=dmean_o3.max(ind);
        
        plot(1,1), hold on % need to have something plotted before turning hold on,
                           % tight_subplot removes ylabels otherwise...

        for dy=1:length(plot_x)
            
            % draw rectangle for each day, with min and max BrO columns as
            % the top and bottom edges
            % position is [x,y,w,h] (bottom left corner)
            
            % plot o3 data first
            if ~isnan(plot_mean_o3)
                rectangle('position',...
                          [plot_x(dy)-box_w/2,plot_min_o3(dy),box_w,plot_max_o3(dy)-plot_min_o3(dy)],...
                          'facecolor','none','edgecolor',[.55 .55 .55])
                rectangle('position',...
                          [plot_x(dy)-box_w/2,plot_mean_o3(dy),box_w,1e12],...
                          'facecolor',[.4 .4 .4],'edgecolor',[.4 .4 .4])
            end
            
            rectangle('position',...
                      [plot_x(dy)-box_w/2,plot_min(dy),box_w,plot_max(dy)-plot_min(dy)],...
                      'facecolor',c_list(2,:),'edgecolor',c_list(2,:))
            rectangle('position',...
                      [plot_x(dy)-box_w/2,plot_mean(dy),box_w,1e12],...
                      'facecolor','k','edgecolor','k')
                                    
        end

        grid on;
        ylim([-0.5,23]*1e13)
        xlim([datenum(yr,3,1),datenum(yr,6,2)])
                
        set(gca, 'XTick', datenum([[yr,3,1];[yr,3,15];[yr,4,1];[yr,4,15];...
                                   [yr,5,1];[yr,5,15];[yr,6,1]]))
        set(gca, 'XTicklabel', {'01/03','15/03','01/04','15/04','01/05','15/05','01/06'})
        
        if yr~=2019
            set(gca,'xticklabel',[])
        else
            xlabel('Date (dd/mm, UTC)')
        end
        
        if yr==2017
            ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
            ylb.Position(1)=ylb.Position(1)-1;
            ylb.Position(2)=ylb.Position(2)-max(dmean_bro.mean);
            
            text(datenum([yr,6,8]),-3e13,'Surface O_3 (ppbv)', 'color','k',...
            'fontsize',fig_fs,'horizontalalignment','center','rotation',90)

        end
        
        text(0.9,0.8,num2str(yr), 'color','k','Units','normalized',...
        'fontsize',fig_fs,'fontweight','bold')

        % create o3 data axes on right side
        text(datenum([yr,6,3]),20e13,'40', 'color','k',...
        'fontsize',fig_fs,'horizontalalignment','left')
        text(datenum([yr,6,3]),10e13,'20', 'color','k',...
        'fontsize',fig_fs,'horizontalalignment','left')
        text(datenum([yr,6,3]),0e13,'0', 'color','k',...
        'fontsize',fig_fs,'horizontalalignment','left')

%         set(gca, 'YMinorTick','on', 'YMinorGrid','on')    
        
    end
    
    
    
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)

    save_pdf(save_figs, 'BrO_o3_dmean')
    
end

%% pca results
if plot_pca

    cutoff=0.25;
    
    eval(['c_list_tmp=' c_scale '(10);'])
    c_list_tmp=c_list_tmp([4,6],:);
    
    figure
    set(gcf, 'Position', [100, 100, 1000, 410]);
    fig_ax = tight_subplot(2,2,[-0.001,0.05],[0.22,0.1],[0.09,0.04]);

    %%% N winds, PC1
    ind=bee_dataset.N_SE_rest==1;
    [vars,coeffs,tot_var]=calc_pca(bee_dataset,ind);

    axes(fig_ax(1))
    barplot_pca(coeffs(:,1),cutoff,tot_var(1),c_list_tmp)
    
    set(gca,'XTickLabel',[])
    ylabel('PC 1 Loading')
    title('{\bfa)} Wind: N','fontweight','normal','FontSize',fig_fs)    
    
    %%% PC2
    axes(fig_ax(3))
    barplot_pca(coeffs(:,2),cutoff,tot_var(2),c_list_tmp)
    
    set(gca,'XTickLabel',vars,'XTickLabelRotation',-45)
    ylabel('PC 2 Loading')
    
    
    %%% SE winds, PC1
    ind=bee_dataset.N_SE_rest==2;
    [~,coeffs,tot_var]=calc_pca(bee_dataset,ind);

    axes(fig_ax(2))
    barplot_pca(coeffs(:,1),cutoff,tot_var(1),c_list_tmp)
    
    set(gca,'XTickLabel',[])
    title('{\bfb)} Wind: SE','fontweight','normal','FontSize',fig_fs)    
    
    %%% SE winds, PC2
    axes(fig_ax(4))
    barplot_pca(coeffs(:,2),cutoff,tot_var(2),c_list_tmp)
    
    set(gca,'XTickLabel',vars,'XTickLabelRotation',-45)
    
           
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)

    save_pdf(save_figs, 'PCA')

    
end
    
    
%% plot relationship to wind speed
if weather_corr
    
    % set up figure
    figure
    set(gcf, 'Position', [100, 100, 1000, 700]);
    fig_ax = tight_subplot(2,3,[0.12,0.05],[0.104,0.06],[0.087,0.04]);
    
    %%%%%% wind speed
%     edges=0:18;
    edges=1;
    
    % northerly winds, mean: 354 deg from gaussian fit, +-30 deg
    axes(fig_ax(1))

    dscatter(bee_dataset.wspd_ms(ind_N), bee_dataset.bro_col(ind_N),'cmap',c_scale)
    hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind_N),bee_dataset.bro_col(ind_N),edges,fig_fs,'N','a',1)
    
    ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
    ylb.Position(1)=-2;
    ylb.Position(2)=bro_lim*1e13/2;
    
    xlabel('Wind speed (m s^{-1})')
    ylim([0,bro_lim]*1e13)
    xlim([0,15.9])   
    
    % southeasterly winds, mean: 123 deg from gaussian fit, +-30 deg  
    axes(fig_ax(2))

    dscatter(bee_dataset.wspd_ms(ind_SE), bee_dataset.bro_col(ind_SE),'cmap',c_scale)
    hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind_SE),bee_dataset.bro_col(ind_SE),edges,fig_fs,'SE','b',1)
    
%     ylabel('BrO part. col. (molec cm^{-2})')
    xlabel('Wind speed (m s^{-1})')
    ylim([0,bro_lim]*1e13)
    xlim([0,15.9])   

    % everything else
    axes(fig_ax(3))

    dscatter(bee_dataset.wspd_ms(ind_rest), bee_dataset.bro_col(ind_rest),'cmap',c_scale)
    hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind_rest),bee_dataset.bro_col(ind_rest),edges,fig_fs,...
                  'other','c',1)
    
%     ylabel('BrO part. col. (molec cm^{-2})')
    xlabel('Wind speed (m s^{-1})')
    ylim([0,bro_lim]*1e13)
    xlim([0,15.9])   
    
    %%%%%% temperature
    edges=1;
    ind_t_rl=~isnan(bee_dataset.T_PWS);
    
    % northerly winds, mean: 354 deg from gaussian fit, +-30 deg
    axes(fig_ax(4))
    ind=(ind_N & ind_t_rl);
    
    dscatter(bee_dataset.T_PWS(ind), bee_dataset.bro_col(ind),'cmap',c_scale)
    hold on, box on
    plot_mean_std(bee_dataset.T_PWS(ind),bee_dataset.bro_col(ind),edges,fig_fs,'N','d',99)
%     plot_fit_line(bee_dataset.T_PWS(ind),bee_dataset.bro_col(ind),fig_fs,'N','d',1)

    ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
    ylb.Position(1)=-53;
    ylb.Position(2)=bro_lim*1e13/2;
    
    xlabel('Temperature (°C)')
    ylim([0,bro_lim]*1e13)
    xlim([-46,6.5])    
    
    % southeasterly winds, mean: 123 deg from gaussian fit, +-30 deg  
    axes(fig_ax(5))
    ind=(ind_SE & ind_t_rl);

    dscatter(bee_dataset.T_PWS(ind), bee_dataset.bro_col(ind),'cmap',c_scale)
    hold on, box on
    plot_mean_std(bee_dataset.T_PWS(ind),bee_dataset.bro_col(ind),edges,fig_fs,'SE','e',99)
%     plot_fit_line(bee_dataset.T_PWS(ind),bee_dataset.bro_col(ind),fig_fs,'SE','e',1)
    
%     ylabel('BrO part. col. (molec cm^{-2})')
    xlabel('Temperature (°C)')
    ylim([0,bro_lim]*1e13)
    xlim([-46,6.5])    

    % everything else
    axes(fig_ax(6))
    ind=(ind_rest & ind_t_rl);

    dscatter(bee_dataset.T_PWS(ind), bee_dataset.bro_col(ind),'cmap',c_scale)
    hold on, box on
    plot_mean_std(bee_dataset.T_PWS(ind),bee_dataset.bro_col(ind),edges,fig_fs,'other','f',99)
%     plot_fit_line(bee_dataset.T_PWS(ind),bee_dataset.bro_col(ind),fig_fs,'other','f',1)
    
%     ylabel('BrO part. col. (molec cm^{-2})')
    xlabel('Temperature (°C)')
    ylim([0,bro_lim]*1e13)
    xlim([-46,6.5])    


    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)

    
    save_pdf(save_figs, 'BrO_weather')
    
end


%% ozone, AOD and aerosols as a function of wind speed
if o3_aer_wspd

    % set up figure
    figure
    set(gcf, 'Position', [100, 100, 1000, 335]);
    fig_ax = tight_subplot(1,3,[0.1,0.08],[0.23,0.08],[0.087,0.04]);
    
    edges=1;
    
    % ozone plot
    axes(fig_ax(1))
    
    % using data paired to BrO
    ind=(~isnan(bee_dataset.wspd_ms_EWS) & ~isnan(bee_dataset.o3_surf));
    
    dscatter(bee_dataset.wspd_ms_EWS(ind),bee_dataset.o3_surf(ind),'cmap',c_scale), hold on, box on
    plot_mean_std(bee_dataset.wspd_ms_EWS(ind),bee_dataset.o3_surf(ind),edges,fig_fs,'','a')

% % %     % using original data
% % %     load('/home/kristof/work/weather_stations/Eureka/EWS_PTU_and_weather_complete.mat')
% % %     ews_data=data;
% % %     load('/home/kristof/work/surface_ozone/surf_o3_hourly_all.mat');
% % %     
% % %     wspd_ews=interp1(ews_data.DateTime,ews_data.WindSpdkmh*10/36,surf_o3_hourly.DateTime,'nearest','extrap');    
% % %     ind=(~isnan(wspd_ews) & ~isnan(surf_o3_hourly.o3_ppb));
% % %     
% % %     dscatter(wspd_ews(ind),surf_o3_hourly.o3_ppb(ind),'cmap',c_scale), hold on, box on
% % %     plot_mean_std(wspd_ews(ind),surf_o3_hourly.o3_ppb(ind),edges,fig_fs,'','a')
    
    xlim([0,15.9])
    ylim([0,45])
    
    xlabel('Wind speed, EWS (m s^{-1})')
    ylb=ylabel('Surface ozone (ppbv)');
    ylb.Position(1)=-2.3;
    
    
    % AOD plot
    axes(fig_ax(2))

    ind=~isnan(bee_dataset.N_SE_rest);
    dscatter(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)),'cmap',c_scale), hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)),edges,fig_fs,'','b')
    xlim([0,15.9])
    ylim(log([0.007,5.1]))

    h = gca;
    set(h,'YTick',log([0.01,0.1,1]))
    set(h,'YTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    h.YAxis.MinorTick = 'on'; 
    h.YAxis.MinorTickValues=log([[1e-3:1e-3:9e-3],[1e-2:1e-2:9e-2],[0.1:0.1:0.9],[1:1:9]]); 
    
    xlabel('Wind speed (m s^{-1})')
    ylabel('AOD')
    
    % aerosol plot
    axes(fig_ax(3))

    ind=(~isnan(bee_dataset.aer_halfmicron) & ~isnan(bee_dataset.N_SE_rest));
    dscatter(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_halfmicron(ind)),'cmap',c_scale), hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_halfmicron(ind)),edges,fig_fs,'','c')
    xlim([0,15.9])
    ylim(log([0.007,5.1]))

    h = gca;
    set(h,'YTick',log([0.01,0.1,1]))
    set(h,'YTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    h.YAxis.MinorTick = 'on'; 
    h.YAxis.MinorTickValues=log([[1e-3:1e-3:9e-3],[1e-2:1e-2:9e-2],[0.1:0.1:0.9],[1:1:9]]); 
    
    xlabel('Wind speed (m s^{-1})')
    ylb=ylabel('d_p > 0.5 \mum (cm^{-3})');
    ylb.Position(1)=-2.6;
    
    
    % figure properties
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
    save_pdf(save_figs, 'o3_aer_wspd')
    
end


%% back trajectories and surface sensitivity

if sens_map
    
    font_correction=0;
    
% % %     % using EWS wind data
% % %     load('/home/kristof/work/weather_stations/Eureka/EWS_PTU_and_weather_complete.mat')
% % %     ews_data=data;
% % %     wspd_ews=interp1(ews_data.DateTime,ews_data.WindSpdkmh*10/36,bee_fp.mean_time,'nearest','extrap');    
    
    % index for later: remove 2015, and remove FP runs where all BrO meas
    % have been filtered out
    fp_good=(bee_fp.mean_time.Year~=2015 & ~isnan(bee_fp.bro_mean_col));

    
    
    % load sensitivity and back trajectory data
    load('/home/kristof/berg/FLEXPART_10.02/BrO_back_runs_v2/BrO_back_runs_v2_3day.mat')
% % %     latitude=1;
% % %     longitude=1;
% % %     sensitivities=1;
% % %     trajectories=1;
% % %     plot_ind=1;

    %%%%%%% condition to select data to plot in top and bottom rows (other than wdir)
    % mean (includes 2015)
%     row1=bee_fp.bro_m==2;
%     row2=bee_fp.bro_m==1;
    % 3 * approx. detection limit
%     row1=bee_fp.bro_mean_col>=1.11e13; 
%     row2=bee_fp.bro_mean_col<1.11e13;
    % by quartile
    qrtiles=prctile(bee_fp.bro_mean_col(fp_good),[25,50,75]);
    row1=bee_fp.bro_mean_col>=qrtiles(3); 
    row2=(bee_fp.bro_mean_col<qrtiles(3) & bee_fp.bro_mean_col>=qrtiles(2));
    row3=bee_fp.bro_mean_col<qrtiles(2);
    
    qrtiles=prctile(bee_dataset.bro_col,[25,50,75]);
    all_row1=bee_dataset.bro_col>=qrtiles(3);
    all_row2=(bee_dataset.bro_col<qrtiles(3) & bee_dataset.bro_col>=qrtiles(2));
    all_row3=bee_dataset.bro_col<qrtiles(2);
    

    %%%%%%%
    
    figure
%     set(gcf, 'Position', [100, 100, 1000, 800]);
%     fig_ax = tight_subplot(2,3,[0.06,0.02],[0.12,0.05],[0.035,0.035]);
    set(gcf, 'Position', [100, 100, 1000, 1000]);
%     fig_ax = tight_subplot(3,3,[0.04,0.02],[0.12,0.05],[0.035,0.035]);
    fig_ax = tight_subplot(3,3,[0.02,0.02],[0.10,0.055],[0.06,0.03]);
    
    % N winds and ?
    axes(fig_ax(1))
    plot_ind=(bee_fp.wdir==1 & row1 & fp_good); 
%     fraction=(sum(bee_dataset.N_SE_rest==1 & all_row1)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;
    
    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'a',fraction);
% %     title('{\bfa)} Wind: N','fontweight','normal','FontSize',fig_fs+font_correction)
    
    % SE winds and ?
    axes(fig_ax(2))
    plot_ind=(bee_fp.wdir==2 & row1 & fp_good); 
%     fraction=(sum(bee_dataset.N_SE_rest==2 & all_row1)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;
    
    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'b',fraction);
% %     title('{\bfb)} Wind: SE','fontweight','normal','FontSize',fig_fs+font_correction)

    % other winds and ?
    axes(fig_ax(3))
    plot_ind=(bee_fp.wdir==3 & row1 & fp_good);
%     fraction=(sum(bee_dataset.N_SE_rest==3 & all_row1)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;

    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'c',fraction);
% %     title('{\bfc)} Wind: other','fontweight','normal','FontSize',fig_fs+font_correction)
    
    %%%%%%%%%%%%%%%%%
    
    % N winds and ?
    axes(fig_ax(4))
    plot_ind=(bee_fp.wdir==1 & row2 & fp_good); 
%     fraction=(sum(bee_dataset.N_SE_rest==1 & all_row2)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;

    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'d',fraction);
% %     title('{\bfd)} Wind: N','fontweight','normal','FontSize',fig_fs+font_correction)
    
    % SE winds and ?
    axes(fig_ax(5))
    plot_ind=(bee_fp.wdir==2 & row2 & fp_good); 
%     fraction=(sum(bee_dataset.N_SE_rest==2 & all_row2)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;

    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'e',fraction);%,...
%                   [2/6,0.08,2/6,0.03]); % set colorbar position
% %     title('{\bfe)} Wind: SE','fontweight','normal','FontSize',fig_fs+font_correction)

    % other winds ?
    axes(fig_ax(6))
    plot_ind=(bee_fp.wdir==3 & row2 & fp_good);
%     fraction=(sum(bee_dataset.N_SE_rest==3 & all_row2)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;

    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'f',fraction);
% %     title('{\bff)} Wind: other','fontweight','normal','FontSize',fig_fs+font_correction)

    %%%%%%%%%%%%%%%%%
    
    % N winds and ?
    axes(fig_ax(7))
    plot_ind=(bee_fp.wdir==1 & row3 & fp_good); 
%     fraction=(sum(bee_dataset.N_SE_rest==1 & all_row3)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;

    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'g',fraction);
% %     title('{\bfg)} Wind: N','fontweight','normal','FontSize',fig_fs+font_correction)
    
    % SE winds and ?
    axes(fig_ax(8))
    plot_ind=(bee_fp.wdir==2 & row3 & fp_good); 
%     fraction=(sum(bee_dataset.N_SE_rest==2 & all_row3)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;

    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'h',fraction,...
                  [2/6+0.015,0.06,2/6,0.02]); % set colorbar position
%                   [2/6,0.08,2/6,0.02]); % set colorbar position
% %     title('{\bfh)} Wind: SE','fontweight','normal','FontSize',fig_fs+font_correction)

    % other winds ?
    axes(fig_ax(9))
    plot_ind=(bee_fp.wdir==3 & row3 & fp_good);
%     fraction=(sum(bee_dataset.N_SE_rest==3 & all_row3)/length(bee_dataset.bro_col))*100;
    fraction=(sum(plot_ind)/sum(fp_good))*100;

    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,'i',fraction);
% %     title('{\bfi)} Wind: other','fontweight','normal','FontSize',fig_fs+font_correction)
    

    %%%tmp
    
    axes(fig_ax(1))
    text(0.5,1.1,'Wind: N', 'color','k','Units','normalized',...
    'fontsize',14,'HorizontalAlignment','center')
    axes(fig_ax(2))
    text(0.5,1.1,'Wind: SE', 'color','k','Units','normalized',...
    'fontsize',14,'HorizontalAlignment','center')
    axes(fig_ax(3))
    text(0.5,1.1,'Wind: other', 'color','k','Units','normalized',...
    'fontsize',14,'HorizontalAlignment','center')
    axes(fig_ax(1))
    text(-0.11,0.5,'BrO > 75^{th} pctl.', 'color','k','Units','normalized',...
    'fontsize',14,'HorizontalAlignment','center','rotation',90)
    axes(fig_ax(4))
    text(-0.11,0.5,'50^{th} pctl. < BrO < 75^{th} pctl.', 'color','k','Units','normalized',...
    'fontsize',14,'HorizontalAlignment','center','rotation',90)
    axes(fig_ax(7))
    text(-0.11,0.5,'BrO < 50^{th} pctl.', 'color','k','Units','normalized',...
    'fontsize',14,'HorizontalAlignment','center','rotation',90)

    
    %%%


    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
%     if save_figs, save_pdf(2, 'sens_map'), end
    if save_figs, error('Resize manually'), end
    % cannot set fig size larger than screen by default -- resize manually
    % to larger than screen, and then 1000x1000 pixels will work
    
end



%% correlation with coarse mode particles and AOD
if plot_ssa
    
    
    figure
%     set(gcf, 'Position', [100, 100, 1000, 370]);
%     fig_ax = tight_subplot(1,3,[0.08,0.05],[0.21,0.113],[0.087,0.04]);
    set(gcf, 'Position', [100, 100, 1000, 700]);
    fig_ax = tight_subplot(2,3,[0.12,0.05],[0.104,0.06],[0.087,0.04]);
        
    ind_ssa=(bee_dataset.aer_halfmicron <= 100);

    plot_type='hm';
    
    switch plot_type
        case 'sm'
            plot_data=bee_dataset.aer_supermicron;
            xlim_end=1;
            x_label='d_p > 1 \mum (cm^{-3})';
        case 'hm'
            plot_data=bee_dataset.aer_halfmicron;
            xlim_end=5.2;
            x_label='d_p > 0.5 \mum (cm^{-3})';
    end
    
    % Northerly winds only
    axes(fig_ax(4))
    ind=(ind_ssa & ind_N);
    dscatter(plot_data(ind), bee_dataset.bro_col(ind),'cmap',c_scale), hold on, box on
    
    plot_fit_line(plot_data(ind),bee_dataset.bro_col(ind),fig_fs,'N','d')    
    
    ylim([0,bro_lim]*1e13)
    xlim([0,xlim_end])    
    xlb=xlabel(x_label);
    xlb.Position(2)=-bro_lim*(9/8)*1e12;
    ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
    ylb.Position(1)=-0.68;
    ylb.Position(2)=bro_lim*1e13/2;
    

    % Southeasterly winds only
    axes(fig_ax(5))
    ind=(ind_ssa & ind_SE);
    dscatter(plot_data(ind), bee_dataset.bro_col(ind),'cmap',c_scale), hold on, box on
    
    plot_fit_line(plot_data(ind),bee_dataset.bro_col(ind),fig_fs,'SE','e')    
    
    ylim([0,bro_lim]*1e13)
    xlim([0,xlim_end])    
    xlb=xlabel(x_label);
    xlb.Position(2)=-bro_lim*(9/8)*1e12;

    % all other wind directions
    axes(fig_ax(6))
    ind=(ind_ssa & ind_rest);
    dscatter(plot_data(ind), bee_dataset.bro_col(ind),'cmap',c_scale), hold on, box on
    
    plot_fit_line(plot_data(ind),bee_dataset.bro_col(ind),fig_fs,'other','f')    
    
    ylim([0,bro_lim]*1e13)
    xlim([0,xlim_end])    
    xlb=xlabel(x_label);
    xlb.Position(2)=-bro_lim*(9/8)*1e12;

    
    %%% AOD corr
    
     
%     set(gca,'XTick',log([0.01,0.1,1]))
%     set(gca,'XMinorTick','on')%log([[0.01:0.01:0.09],[0.1:0.1:0.9],[1:1:2]]))
%     tmp=gca;
%     tmp.XAxis.MinorTickValues = log([[0.007,0.008,0.009],[0.01:0.01:0.09],[0.1:0.1:0.9],[1:1:7]]);
%     set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    axes(fig_ax(1))
    dscatter(log(bee_dataset.aer_ext(ind_N)),bee_dataset.bro_col(ind_N),'cmap',c_scale);hold on, box on
    plot_mean_std(log(bee_dataset.aer_ext(ind_N)),bee_dataset.bro_col(ind_N),1,fig_fs,'N','a')
    xlb=xlabel('AOD');
%     xlb.Position(2)=-11e12;
    ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
    ylb.Position(1)=-5.7;
    ylb.Position(2)=bro_lim*1e13/2;
    
    ylim([0,bro_lim]*1e13)
    xlim([-4.8,2])

    set(gca,'XTick',log([0.01,0.1,1]))
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    axes(fig_ax(2))
    dscatter(log(bee_dataset.aer_ext(ind_SE)),bee_dataset.bro_col(ind_SE),'cmap',c_scale);hold on, box on
    plot_mean_std(log(bee_dataset.aer_ext(ind_SE)),bee_dataset.bro_col(ind_SE),1,fig_fs,'SE','b')
    xlb=xlabel('AOD');
%     xlb.Position(2)=-11e12;
    ylim([0,bro_lim]*1e13)
    xlim([-4.8,2])

    set(gca,'XTick',log([0.01,0.1,1]))
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    axes(fig_ax(3))
    dscatter(log(bee_dataset.aer_ext(ind_rest)),bee_dataset.bro_col(ind_rest),'cmap',c_scale);hold on, box on
    plot_mean_std(log(bee_dataset.aer_ext(ind_rest)),bee_dataset.bro_col(ind_rest),1,fig_fs,'other','c')
    xlb=xlabel('AOD');
%     xlb.Position(2)=-11e12;
    ylim([0,bro_lim]*1e13)
    xlim([-4.8,2])
   
    set(gca,'XTick',log([0.01,0.1,1]))
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    
    
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
    
    save_pdf(save_figs, 'aer_corr')
    
    
end


%% BrO columns as a function of SI contact 
if si_contact_log
    
    % 1: plot interpolated contact value for individual BrO columns
    % 0: plot contact value for each FP run as a function of the mean BrO column
    plot_interp=1; 
    
    % manually position x label multiplier to get it out of the way
    labelmult_x=8.9*1e13;
    labelmult_y=-0.4*1e13;

    if plot_log_si

        if strcmp(btraj_len,'3')
            lim_si=[50,9.4e5];
        elseif strcmp(btraj_len,'4')
            lim_si=[4e2,1.4e6];
        elseif strcmp(btraj_len,'5')
            lim_si=[2e3,1.5e6];
        end
        
        si_x_lim=log(lim_si);
        edges=1;
        flip_ind=3;
        
        xtick_pos=log([10,1e2,1e3,1e4,1e5,1e6]);
        xtickm_pos=log([[10:10:90],[1e2:1e2:9e2],[1e3:1e3:9e3],...
                        [1e4:1e4:9e4],[1e5:1e5:9e5],[1e6:1e6:9e6]]);
        xtick_label={'10','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'};
        
    else
        si_x_lim=[10,8.1e5];
        edges=[0:13]*1e13;   
        flip_ind=0;
    end
%     v_edges=[0:0.5:9]*1e13;   
    
    %%%%%%%% FYSI %%%%%%%%%    
    if plot_interp % plot interpolated contact value for individual BrO columns
        
        plot_y=bee_dataset.bro_col;

        if plot_log_si
            % must use log of data to get the point desity right (in dscatter)
            eval(['plot_x=log(bee_dataset.FYSI_' btraj_len 'day);']);
        else
            eval(['plot_x=bee_dataset.FYSI_' btraj_len 'day;']);
        end
        
        ind_fp=(~isnan(plot_x) & plot_x>0 & ~isnan(plot_y)); %skip few zero values
        ind_N_tmp=(ind_fp & ind_N);
        ind_SE_tmp=(ind_fp & ind_SE);
        ind_other_tmp=(ind_fp & ind_rest);
        
    else % plot contact value for each FP run, as a function of the mean BrO column
        
        plot_y=bee_fp.bro_mean_col;
        if plot_log_si
            plot_x=log(bee_fp.FYSI_3day);
        else
            plot_x=bee_fp.FYSI_3day;
        end
        
        ind_fp=(~isnan(plot_x) & plot_x>0 & ~isnan(plot_y)); %skip few zero values
        ind_N_tmp=(ind_fp & bee_fp.wdir==1);
        ind_SE_tmp=(ind_fp & bee_fp.wdir==2);
        ind_other_tmp=(ind_fp & bee_fp.wdir>=3);
        
    end
        

    figure
    set(gcf, 'Position', [100, 100, 1000, 700]);
%     fig_ax = tight_subplot(2,3,[0.12,0.05],[0.114,0.06],[0.087,0.04]);
    fig_ax = tight_subplot(2,3,[0.12,0.05],[0.104,0.06],[0.087,0.04]);

    % Northerly winds only
    axes(fig_ax(1))
    dscatter(plot_x(ind_N_tmp), plot_y(ind_N_tmp),'cmap',c_scale), hold on, box on

    plot_mean_std(plot_x(ind_N_tmp),plot_y(ind_N_tmp),edges,fig_fs,'N','a',flip_ind)
%     plot_vertical_mean_std(plot_x(ind_N_tmp),plot_y(ind_N_tmp),v_edges)
%     plot_fit_line(plot_x(ind_N_tmp),plot_y(ind_N_tmp),fig_fs,'N','a')

    ylim([0,bro_lim]*1e13)
    xlim(si_x_lim)    
    ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
    ylb.Position(1)=log(lim_si(1))-0.1313315*log(lim_si(2)/lim_si(1));
    ylb.Position(2)=bro_lim*1e13/2;
    
    xlb=xlabel('FYI sensitivity (s)');
    xlb.Position(2)=-bro_lim*(11/8)*1e12;

    if plot_log_si
        h = gca;
        set(h,'XTick',xtick_pos)
        set(h,'XTickLabel',xtick_label)
        h.XAxis.MinorTick = 'on'; 
        h.XAxis.MinorTickValues=xtickm_pos; 
    else
%         set(gca,'Xtick',[0,2,4,6,8]*1e13)
%         set(gca,'XtickLabel',[0,2,4,6,8])
%         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
    end    
    
    % Southeasterly winds only
    axes(fig_ax(2))
    dscatter(plot_x(ind_SE_tmp), plot_y(ind_SE_tmp),'cmap',c_scale), hold on, box on

    plot_mean_std(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),edges,fig_fs,'SE','b',flip_ind)
%     plot_vertical_mean_std(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),v_edges)
%     plot_fit_line(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),fig_fs,'SE','b')

    ylim([0,bro_lim]*1e13)
    xlim(si_x_lim)    
    xlb=xlabel('FYI sensitivity (s)');
    xlb.Position(2)=-bro_lim*(11/8)*1e12;
    
    if plot_log_si
        h = gca;
        set(h,'XTick',xtick_pos)
        set(h,'XTickLabel',xtick_label)
        h.XAxis.MinorTick = 'on'; 
        h.XAxis.MinorTickValues=xtickm_pos; 
    else
%         set(gca,'Xtick',[0,2,4,6,8]*1e13)
%         set(gca,'XtickLabel',[0,2,4,6,8])
%         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
    end    

    
    % all other wind directions
    axes(fig_ax(3))
    dscatter(plot_x(ind_other_tmp), plot_y(ind_other_tmp),'cmap',c_scale), hold on, box on

    plot_mean_std(plot_x(ind_other_tmp),plot_y(ind_other_tmp),edges,fig_fs,'other','c',flip_ind)
%     plot_vertical_mean_std(plot_x(ind_other_tmp),plot_y(ind_other_tmp),v_edges)
%     plot_fit_line(plot_x(ind_other_tmp),plot_y(ind_other_tmp),fig_fs,'other','c')

    ylim([0,bro_lim]*1e13)
    xlim(si_x_lim)    
    xlb=xlabel('FYI sensitivity (s)');
    xlb.Position(2)=-bro_lim*(11/8)*1e12;
    
    if plot_log_si
        h = gca;
        set(h,'XTick',xtick_pos)
        set(h,'XTickLabel',xtick_label)
        h.XAxis.MinorTick = 'on'; 
        h.XAxis.MinorTickValues=xtickm_pos; 
    else
%         set(gca,'Xtick',[0,2,4,6,8]*1e13)
%         set(gca,'XtickLabel',[0,2,4,6,8])
%         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
    end    
        
    
    %%%%%%%% All SI %%%%%%%%%    
    if plot_interp % plot interpolated contact value for individual BrO columns
        
        plot_y=bee_dataset.bro_col;

        if plot_log_si
            % must use log of data to get the point desity right (in dscatter)
            eval(['plot_x=log(bee_dataset.FYSI_' btraj_len 'day',...
                  ' + bee_dataset.MYSI_' btraj_len 'day);']);
%             si_x_lim=log([3e9,1.1e14]);
        else
            eval(['plot_x=bee_dataset.FYSI_' btraj_len 'day',...
                  ' + bee_dataset.MYSI_' btraj_len 'day;']);
        end
        
%         ind_fp=(~isnan(plot_x) & ~isnan(plot_y));
        ind_fp=(~isnan(plot_x) & bee_dataset.FYSI_3day>0 & ~isnan(plot_y)); %skip 0s for 3 day
        ind_N_tmp=(ind_fp & ind_N);
        ind_SE_tmp=(ind_fp & ind_SE);
        ind_other_tmp=(ind_fp & ind_rest);
        
    else % plot contact value for each FP run, as a function of the mean BrO column
        
        plot_y=bee_fp.bro_mean_col;
        plot_x=bee_fp.FYSI_3day;
        
        ind_N_tmp=(bee_fp.wdir==1);
        ind_SE_tmp=(bee_fp.wdir==2);
        ind_other_tmp=(bee_fp.wdir>=3);
        
    end
        

    % Northerly winds only
    axes(fig_ax(4))
    dscatter(plot_x(ind_N_tmp), plot_y(ind_N_tmp),'cmap',c_scale), hold on, box on
    
%     ratio=bee_dataset.MYSI_3day./(bee_dataset.FYSI_3day+bee_dataset.MYSI_3day);    
%     plot(plot_x(ind_N_tmp & ratio>0.9 & bee_dataset.o3_surf<15), plot_y(ind_N_tmp & ratio>0.9 & bee_dataset.o3_surf<15),'y.')

    plot_mean_std(plot_x(ind_N_tmp),plot_y(ind_N_tmp),edges,fig_fs,'N','d',flip_ind)
%     plot_vertical_mean_std(plot_x(ind_N_tmp),plot_y(ind_N_tmp),v_edges)
%     plot_fit_line(plot_x(ind_N_tmp),plot_y(ind_N_tmp),fig_fs,'N','a')
    
    ylim([0,bro_lim]*1e13)
    xlim(si_x_lim)    
    ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
    ylb.Position(1)=log(lim_si(1))-0.1313315*log(lim_si(2)/lim_si(1));
    ylb.Position(2)=bro_lim*1e13/2;
    
    xlb=xlabel('Total ice sensitivity (s)');
    xlb.Position(2)=-bro_lim*(11/8)*1e12;

    if plot_log_si
        h = gca;
        set(h,'XTick',xtick_pos)
        set(h,'XTickLabel',xtick_label)
        h.XAxis.MinorTick = 'on'; 
        h.XAxis.MinorTickValues=xtickm_pos; 
    else
%         set(gca,'Xtick',[0,2,4,6,8]*1e13)
%         set(gca,'XtickLabel',[0,2,4,6,8])
%         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
    end    
    
    % Southeasterly winds only
    axes(fig_ax(5))
    dscatter(plot_x(ind_SE_tmp), plot_y(ind_SE_tmp),'cmap',c_scale), hold on, box on

    plot_mean_std(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),edges,fig_fs,'SE','e',flip_ind)
%     plot_vertical_mean_std(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),v_edges)
%     plot_fit_line(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),fig_fs,'SE','b')

    ylim([0,bro_lim]*1e13)
    xlim(si_x_lim)    
    xlb=xlabel('Total ice sensitivity (s)');
    xlb.Position(2)=-bro_lim*(11/8)*1e12;
    
    if plot_log_si
        h = gca;
        set(h,'XTick',xtick_pos)
        set(h,'XTickLabel',xtick_label)
        h.XAxis.MinorTick = 'on'; 
        h.XAxis.MinorTickValues=xtickm_pos; 
    else
%         set(gca,'Xtick',[0,2,4,6,8]*1e13)
%         set(gca,'XtickLabel',[0,2,4,6,8])
%         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
    end    

    
    % all other wind directions
    axes(fig_ax(6))
    dscatter(plot_x(ind_other_tmp), plot_y(ind_other_tmp),'cmap',c_scale), hold on, box on

    plot_mean_std(plot_x(ind_other_tmp),plot_y(ind_other_tmp),edges,fig_fs,'other','f',flip_ind)
%     plot_vertical_mean_std(plot_x(ind_other_tmp),plot_y(ind_other_tmp),v_edges)
%     plot_fit_line(plot_x(ind_other_tmp),plot_y(ind_other_tmp),fig_fs,'other','c')

    ylim([0,bro_lim]*1e13)
    xlim(si_x_lim)    
    xlb=xlabel('Total ice sensitivity (s)');
    xlb.Position(2)=-bro_lim*(11/8)*1e12;
    
    if plot_log_si
        h = gca;
        set(h,'XTick',xtick_pos)
        set(h,'XTickLabel',xtick_label)
        h.XAxis.MinorTick = 'on'; 
        h.XAxis.MinorTickValues=xtickm_pos; 
    else
%         set(gca,'Xtick',[0,2,4,6,8]*1e13)
%         set(gca,'XtickLabel',[0,2,4,6,8])
%         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
    end     
     
    
% %     %%%%%%%% MYSI %%%%%%%%%
% %     if plot_interp % plot using all BrO cols and fix Si contact claculated for each
% %         
% %         plot_y=bee_dataset.bro_col;
% %         
% %         if plot_log_si
% %             % must use log of data to get the point desity right (in dscatter)
% %             eval(['bee_dataset.MYSI_' btraj_len 'day(bee_dataset.MYSI_' btraj_len 'day<90)=90;']);
% %             eval(['plot_x=log(bee_dataset.MYSI_' btraj_len 'day);']);
% %         else
% %             eval(['plot_x=bee_dataset.MYSI_' btraj_len 'day;']);
% %         end
% % 
% % %         ind_fp=(~isnan(bee_dataset.MYSI_3day) & bee_dataset.MYSI_3day>0);
% %         ind_fp=~isnan(plot_x);
% %         
% %         ind_N_tmp=(ind_fp & ind_N);
% %         ind_SE_tmp=(ind_fp & ind_SE);
% %         ind_other_tmp=(ind_fp & ind_rest);
% %         
% %     else % plot using FP runs, with calculated mean BrO for each
% %         
% %         plot_y=bee_fp.bro_mean_col;
% %         plot_x=bee_fp.MYSI_3day;
% %         
% %         ind_N_tmp=(bee_fp.wdir==1);
% %         ind_SE_tmp=(bee_fp.wdir==2);
% %         ind_other_tmp=(bee_fp.wdir>=3);
% %         
% %     end
% %     
% %     if plot_log_si
% %         si_x_lim2=log([90,1.02e14]);
% %         edges=1;
% %     else
% %         si_x_lim2=[2e7,1.02e14];
% %         edges=[0:0.5:5]*1e13;
% %     end
% %     
% %     % Northerly winds only
% %     axes(fig_ax(4))
% %     dscatter(plot_x(ind_N_tmp), plot_y(ind_N_tmp),'cmap',c_scale), hold on, box on
% % 
% %     plot_mean_std(plot_x(ind_N_tmp),plot_y(ind_N_tmp),edges,fig_fs,'N','d',flip_ind)
% % %     plot_vertical_mean_std(plot_x(ind_N_tmp),plot_y(ind_N_tmp),v_edges)
% % %     plot_fit_line(plot_x(ind_N_tmp),plot_y(ind_N_tmp),fig_fs,'N','d')
% % 
% %     ylim([0,bro_lim]*1e13)
% %     xlim(si_x_lim2)    
% %     xlb=xlabel('MYSI sensitivity (s)');
% %     xlb.Position(2)=-1e13;
% %     ylb=ylabel('BrO part. col. (molec cm^{-2})'); 
% %     ylb.Position(1)=log(2.3);
% %     ylb.Position(2)=bro_lim*1e13/2;
% %     
% %     if plot_log_si
% %         set(gca,'XTick',log([1e2,1e6,1e10,1e14]))
% %         set(gca,'XTickLabel',{'10^{2}','10^{6}','10^{10}','10^{14}'})
% %     else
% %         set(gca,'Xtick',[0,2,4,6,8]*1e13)
% %         set(gca,'XtickLabel',[0,2,4,6,8])
% %         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
% %     end    
% % 
% %     
% %     % Southeasterly winds only
% %     axes(fig_ax(5))
% %     dscatter(plot_x(ind_SE_tmp), plot_y(ind_SE_tmp),'cmap',c_scale), hold on, box on
% % 
% %     plot_mean_std(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),edges,fig_fs,'SE','e',flip_ind)
% % %     plot_vertical_mean_std(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),v_edges)
% % %     plot_fit_line(plot_x(ind_SE_tmp),plot_y(ind_SE_tmp),fig_fs,'SE','e')
% % 
% %     ylim([0,bro_lim]*1e13)
% %     xlim(si_x_lim2)    
% %     xlb=xlabel('MYSI sensitivity (s)');
% %     xlb.Position(2)=-1e13;
% %     
% %     if plot_log_si
% %         set(gca,'XTick',log([1e2,1e6,1e10,1e14]))
% %         set(gca,'XTickLabel',{'10^{2}','10^{6}','10^{10}','10^{14}'})
% %     else
% %         set(gca,'Xtick',[0,2,4,6,8]*1e13)
% %         set(gca,'XtickLabel',[0,2,4,6,8])
% %         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
% %     end    
% % 
% %     
% %     % all other wind directions
% %     axes(fig_ax(6))
% %     dscatter(plot_x(ind_other_tmp), plot_y(ind_other_tmp),'cmap',c_scale), hold on, box on
% % 
% %     plot_mean_std(plot_x(ind_other_tmp),plot_y(ind_other_tmp),edges,fig_fs,'other','f',flip_ind)
% % %     plot_vertical_mean_std(plot_x(ind_other_tmp),plot_y(ind_other_tmp),v_edges)
% % %     plot_fit_line(plot_x(ind_other_tmp),plot_y(ind_other_tmp),fig_fs,'other','f')
% % 
% %     ylim([0,bro_lim]*1e13)
% %     xlim(si_x_lim2)    
% %     xlb=xlabel('MYSI sensitivity (s)');
% %     xlb.Position(2)=-1e13;
% %     
% %     if plot_log_si
% %         set(gca,'XTick',log([1e2,1e6,1e10,1e14]))
% %         set(gca,'XTickLabel',{'10^{2}','10^{6}','10^{10}','10^{14}'})
% %     else
% %         set(gca,'Xtick',[0,2,4,6,8]*1e13)
% %         set(gca,'XtickLabel',[0,2,4,6,8])
% %         text(labelmult_x,labelmult_y,'\times10^{13}', 'color','k');
% %     end    
    
      
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
    
    if plot_log_si
        save_pdf(save_figs, 'SI_sens_log')
    else
        save_pdf(save_figs, 'SI_sens')
    end
    
    
    
   
    
end

if presentation_plots
    
    figure 
    set(gcf, 'Position', [100, 100, 360, 330]);
    fig_ax = tight_subplot(1,1,[0.07,0.07],[0.2,0.06],[0.22,0.08]);
    
    axes(fig_ax(1))
    
    ind=~isnan(bee_dataset.N_SE_rest);
    dscatter(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)),'cmap',c_scale), hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)),1)
    xlim([0,15])
    ylim(log([0.01,5]))

    set(gca,'YTick',log([0.01,0.1,1]))
    set(gca,'YTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    xlabel('Wind speed (m s^{-1})')
    ylabel('AOD')
    
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
    save_pdf(save_figs, 'presentations/aod_wspd')
    
    %%%%%%%%%%%%%%%%%
    
    figure 
    set(gcf, 'Position', [100, 100, 360, 330]);
    fig_ax = tight_subplot(1,1,[0.07,0.07],[0.2,0.06],[0.22,0.08]);
    
    axes(fig_ax(1))
    
    ind=(~isnan(bee_dataset.aer_halfmicron) & ~isnan(bee_dataset.N_SE_rest));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind),'cmap',c_scale), hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind),1)
    xlim([0,15])
    ylim([0,5])

    xlabel('Wind speed (m s^{-1})')
    ylabel('d_p > 0.5 \mum (cm^{-3})')

    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
    save_pdf(save_figs, 'presentations/ssa_wspd')
  
    %%%%%%%%%%%%%%%%%
    
    figure 
    set(gcf, 'Position', [100, 100, 360, 330]);
    fig_ax = tight_subplot(1,1,[0.07,0.07],[0.2,0.06],[0.22,0.08]);
    
    axes(fig_ax(1))
    
    ind=(~isnan(bee_dataset.aer_halfmicron) & ~isnan(bee_dataset.N_SE_rest));
    dscatter(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_halfmicron(ind)),'cmap',c_scale), hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_halfmicron(ind)),1)
    xlim([0,15])
    ylim(log([0,5]))

    set(gca,'YTick',log([0.01,0.1,1]))
    set(gca,'YTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    xlabel('Wind speed (m s^{-1})')
    ylabel('d_p > 0.5 \mum (cm^{-3})')

    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
    save_pdf(save_figs, 'presentations/ssa_wspd_log')
    
end

end

function save_pdf(save_figs, fname)

    if save_figs
        
        h=gcf;

        set(h,'Units','Inches');

        pos = get(h,'Position');

        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

        if save_figs==1
            % pdf images
            f_out=['/home/kristof/work/documents/paper_bro/figures/' fname '.pdf'];
            print(h,f_out,'-dpdf','-r0')
        
        elseif save_figs==2
            % jpg images
            f_out=['/home/kristof/work/documents/paper_bro/figures/' fname '.jpg'];
            print(h,f_out,'-djpeg','-r200','-opengl')
%             % png images
%             f_out=['/home/kristof/work/documents/paper_bro/figures/' fname '.png'];
%             print(h,f_out,'-dpng','-r400','-opengl')
        elseif save_figs==3
            % bitmap pdf images
            %%% sens_map size: 
            %%% set(gcf, 'Position', [100, 100, 1000, 1000]);
            %%% AFTER position is set: copy-paste first 4 lines (h=gcf, etc)
            %%% f_out=['/home/kristof/work/documents/paper_bro/figures/sens_map.pdf'];
            f_out=['/home/kristof/work/documents/paper_bro/figures/' fname '.pdf'];
            print(h,f_out,'-dpdf','-r200','-opengl')
            
        end
        
        pause(2)
        close(h)
        
    else
        return
    end

end


function ll=plot_yearly_box(data,box_group,a,nbars,c_list,box_lw,box_outlier,do_legend)

    % get position of each bar
    box_pos=sort([1:a:1+a*(nbars-1),2:a:2+a*(nbars-1),3:a:3+a*(nbars-1)]);

    % set tick position as the meanof each group
    tick_pos=[];
    for i=1:3
        tick_pos=[tick_pos,mean(box_pos((nbars)*(i-1)+1:(nbars)*i))];
    end

    % do box plot
    boxplot(data,box_group,'positions',box_pos,'colors',c_list,...
            'jitter',0.5,'symbol','.')
    hold on, box on
        
    % plot mean as well
    group_ind=unique(box_group);
    c_ind=repmat(1:nbars,1,length(group_ind)/nbars);
    for i=1:length(group_ind)
        plot(box_pos(i), nanmean(data(box_group==group_ind(i))),'x','color',c_list(c_ind(i),:))
    end
        
    % some formatting
    set(findobj(gca,'type','line'),'linew',box_lw)
    set(findall(gca,'tag','Outliers'),'MarkerSize',box_outlier);
    
    set(gca,'xtick',tick_pos)
    set(gca,'xticklabel',{'March','April','May'})
    
    xlim([box_pos(1)-0.15,box_pos(end)+0.15])
    
    set(gca, 'YGrid', 'on')

    % legend for one figure only -- position set later
    if do_legend
        ll=legend(flipud(findall(gca,'Tag','Box')),...
               {'2016','2017','2018','2019'},...
               'Orientation','horizontal','location','northeast');
    end
end

function plot_mean_std(xx,yy,edges,text_size,wdir_str,subplot_id,fliplabels)

    % select where to put labels
    y_wnd=0.92; % top corners
    y_let=0.92; 

    if nargin<7 || fliplabels==0 % wind on left, letter on right
        
        x1=0.05;
        x1_align='left';
        x2=0.95;
        x2_align='right';
        
    elseif fliplabels==1 % wind on right, letter on left

        x1=0.95;
        x1_align='right';
        x2=0.05;
        x2_align='left';
        
    elseif fliplabels==2 % same as 0, move wind slightly to the right

        x1=0.15;
        x1_align='left';
        x2=0.95;
        x2_align='right';

    elseif fliplabels==3 % wind and letter both in top left

        x1=0.05;
        x1_align='left';
        y_wnd=0.8;
        x2=0.05;
        x2_align='left';
        
    elseif fliplabels==99
        
        [~, ~, R2] = line_fit(xx,yy);
        
        x1=0.95;
        x1_align='right';
        x2=0.05;
        x2_align='left';
        
        text(x1,0.8,sprintf('{R^2}=%.2f', R2),'color','k','Units','normalized',...
            'fontsize',text_size,'HorizontalAlignment',x1_align)
        
    end

    if length(edges)==1
        % use deciles
        edges=prctile(xx,0:10:100);
    end
    
    % get mean of data
    tmp_mean=NaN(1,length(edges)-1);
    tmp_std=NaN(1,length(edges)-1);
    
    for i=1:length(edges)-1
        
        data=yy(xx>=edges(i) & xx < edges(i+1));
        
        if length(data)>=10
            tmp_mean(i)=nanmean(data);
            tmp_std(i)=nanstd(data);
        end
    end    
    
    plot_x=edges(1:end-1)+diff(edges)/2;
    
    % plot mean and std
    plot(plot_x,tmp_mean,'k-','linewidth',2.5)
    plot(plot_x,tmp_mean+tmp_std,'k:','linewidth',2.5)
    plot(plot_x,tmp_mean-tmp_std,'k:','linewidth',2.5)   

    if nargin>3
        % add text on plot
        if ~isempty(wdir_str)
            text(x1,y_wnd,sprintf(['Wind: ' wdir_str]), 'color','k','Units','normalized',...
                'fontsize',text_size,'HorizontalAlignment',x1_align)
        end
        
        text(x2,y_let,['{\bf' subplot_id ')}'], 'color','k','Units','normalized',...
            'fontsize',text_size,'HorizontalAlignment',x2_align)
    end
%     tmp=corrcoef(exp(xx),yy);
%     tmp=corrcoef(xx,yy);
%     disp(tmp(1,2)^2)
    
end

function plot_vertical_mean_std(xx,yy,edges)

    % get mean of data
    tmp_mean=[];
    tmp_std=[];
    
    for i=1:length(edges)-1
        tmp_mean=[tmp_mean,nanmean(xx(yy>=edges(i) & yy < edges(i+1)))];
        tmp_std=[tmp_std,nanstd(xx(yy>=edges(i) & yy < edges(i+1)))];
    end    
    
    plot_y=edges(1:end-1)+diff(edges)/2;
    
    % plot mean and std
    plot(tmp_mean,plot_y,'-','color',[.7 .7 .7],'linewidth',2)
    plot(tmp_mean+tmp_std,plot_y,'--','color',[.7 .7 .7],'linewidth',2)
    plot(tmp_mean-tmp_std,plot_y,'--','color',[.7 .7 .7],'linewidth',2)   

end

function plot_fit_line(xx,yy,text_size,wdir_str,subplot_id,fliplabels)

    % select where to put labels
    if nargin==5
        
        x1=0.05;
        x1_align='left';
        x2=0.95;
        x2_align='right';
        
    elseif fliplabels

        x1=0.95;
        x1_align='right';
        x2=0.05;
        x2_align='left';
        
    end

    % get line fit
    [slope, y_int, R2] = line_fit(xx,yy);
    
    % plot the results
    plot_x=[min(xx)-max(abs(xx))*0.04:max(xx)/20:max(xx)+max(abs(xx))*0.05];
    plot(plot_x,slope(1)*plot_x+y_int(1),'k','linewidth',2)

    % add text
    text(x1,0.92,sprintf(['Wind: ' wdir_str]), 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment',x1_align)

    text(x1,0.8,sprintf('{R^2}=%.2f', R2),'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment',x1_align)
    
    text(x2,0.92,['{\bf' subplot_id ')}'], 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment',x2_align)
        
    
end

function plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind,...
                       subplot_id,fraction,cbar_pos)

    if nargin==7, cbar_pos=[]; end

    if length(plot_ind)~=size(sensitivities,3)
        error('plot_ind must be a logical, same length as total number of FP runs')
    end
    
    load coast;
    ax = worldmap([61,90], [-180,180]);
    geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])
    hold on

    % subplot lalel
    text(0.92,0.94,['{\bf' subplot_id ')}'], 'color','k','Units','normalized',...
        'fontsize',14,'HorizontalAlignment','left')
    
    % percentage of all data in given subplot
    text(0.92,0.1,sprintf('%0.1f%%', fraction), 'color','k','Units','normalized',...
        'fontsize',14,'HorizontalAlignment','center','Rotation',45)
    
    sensitivities(sensitivities==0)=NaN;

    sens_tmp=nanmean(sensitivities(:,:,plot_ind),3);

    % replace smallest/largest values
    max_val_lin=1300;

    sens_tmp(sens_tmp<0.1)=0.1;
    sens_tmp(sens_tmp>max_val_lin)=max_val_lin;

    % plot on log scale
    sens_tmp=log10(sens_tmp);

    %%% only works if plots are separate; need common colorbar here
% %     % replace max log value with next largest integer if difference is <0.3
% %     [max_val,max_ind] = max(sens_tmp(:));
% %     if max_val-floor(max_val) >0.7
% %         sens_tmp(max_ind)=ceil(max_val);
% %         max_val=ceil(max_val);
% %     end
    max_val=log10(max_val_lin);

    % plot sensitivity
    surfm(latitude,longitude,sens_tmp,'facecolor', 'interp','facealpha',0.75);
    colormap(flipud(hot))

    if ~isempty(cbar_pos)
% % %         surfm([50:10:80],[50:10:80],magic(4))
% % %         max_val=16;
        % add colorbar with manual labels (convert back to lin space)
        cb=colorbar('south','position',cbar_pos);
        cb_lim=-1:ceil(max_val);
        set(cb,'YTick',cb_lim)
        set(cb,'YTick',cb_lim,'YTickLabel',cellstr(num2str(power(10,cb_lim)')))
        ylabel(cb,'Mean surface sensitivity (s)','FontSize',12)
    end

    setm(gca,'MLineLocation',30,'PLineLocation',10,'MLabelLocation',90,...
         'PLabelLocation',10,'MLabelParallel',65,'PLabelMeridian',0,...
         'FontSize',10)

    to_plot=find(plot_ind==1); 
     
    for i=to_plot'
        data=trajectories(trajectories.index==i,:);
        plotm(data.lat,data.lon,'c','linewidth',0.6)
    end

    plotm(80.053, -86.416, 'kp','markerfacecolor','k','markersize',12)

end

function [vars,coeffs,tot_var]=calc_pca(bee_dataset,ind)

%     vars={'BrO','AOD','Aer','V','T','\DeltaT_{200m}','\DeltaT_{610m}','PBLH','P','\DeltaP',...
%           'FYI','MYI'};
%%%   add PBLH if necessary -- calculations are not great for new GRAW
%%%   sondes (2019 data)
%%%   can add pressure tendency as well -- not correlated with anything

    vars={'BrO','AOD','Aer','V','T','\DeltaT_{200m}','\DeltaT_{610m}','P',...
          'FYI','MYI'};

    % create data arrays for pca
    data_in=[bee_dataset.bro_col(ind),...
             bee_dataset.aer_ext(ind),...
             bee_dataset.aer_halfmicron(ind),...
             bee_dataset.wspd_ms(ind),...
             bee_dataset.T_PWS(ind),...
             bee_dataset.sonde_dT(ind),...
             bee_dataset.sonde_T_200(ind)-bee_dataset.sonde_T_0(ind),...
             bee_dataset.P_PWS(ind),...
             bee_dataset.FYSI_3day(ind),...
             bee_dataset.MYSI_3day(ind),...
             ];

    % do PCA

    % center (mean=0) and standardize (std=1) (setting 'VariableWeights','variance'
    % in pca doesn't actually standardize the data for some reason)
    data_mean=nanmean(data_in);
    data_std=nanstd(data_in);

    for i=1:size(data_in,2)
        data_in(:,i)=data_in(:,i)-data_mean(i);
        data_in(:,i)=data_in(:,i)/data_std(i);
    end

    % calculate pca coefficients (loadings), and get the principal component
    % variances and percentage of the total variance explained by each component
    %%% each column of coeffs corresponds to one principal component
    [coeffs,~,~,~,tot_var] = pca(data_in,'Centered',false);


end

function barplot_pca(coeffs,cutoff,tot_var,c_list)

        % loading greater than cutoff, for coloring bars
        ind=abs(coeffs)>cutoff;

        % flip so BrO is consistently positive (easier to interpret
        % results -- sign is arbitrary anyway)
        if coeffs(1)<0, coeffs=-coeffs; end

        % plot loadings above and below the cutoff separately
        x=1:size(coeffs,1);

        y1=coeffs;
        y1(~ind)=0; % > cutoff
        y2=coeffs;
        y2(ind)=0; % < cutoff

        b1=bar(x,y1,'FaceColor',c_list(1,:),'Edgecolor','k'); hold on
        b2=bar(x,y2,'FaceColor',c_list(2,:),'Edgecolor','k');

%         b1.FaceAlpha = 0.5;

        ylim([-0.7,0.7])
        xlim([0.5,size(coeffs,1)+0.5])

        % add percent of variance explained to plot
        text(0.97,0.15,[num2str(round(tot_var,1)) '%'],...
            'color','k','Units','normalized','HorizontalAlignment','right')

end
