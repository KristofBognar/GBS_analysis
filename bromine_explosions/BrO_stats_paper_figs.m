function BrO_stats_paper_figs()
%BRO_STATS_PAPER_FIGS to plot figures for paper describing BrO statistics
%in Eureka
%
% This code contains modified versions of the code from several other
% functions -- collected here to easily update/reproduce figures in the
% paper.
%
% Some data used in plotting are not the usual output files, but copies
% saved in /home/kristof/work/documents/paper_bro (to track data versions)


%% setup

% what to plot
windrose=0;
plot_box=0;
bro_dailymean=0;
weather_corr=0;
sens_meanBrO=0;
plot_ssa=1;
si_contact=0;

% 0 to not save, 1 to save as pdf, 2 to save as jpg
save_figs=0;

% uniform look
fig_fs=14; % font size on figures
fig_font='Arial'; % font for plotting

% line with for box plot, and markers
box_lw=1.5;
box_outlier=4; % outlier size for box plots
    

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
     
m_names={'March','April','May'};
letters='abcdefghijklmnopqrstuvwxyz';


%% load/filter BEE dataset
load('/home/kristof/work/documents/paper_bro/data/BEE_dataset_all.mat')

% exclude 2015 (intermittent measurements)
bee_dataset(bee_dataset.times.Year==2015,:)=[];

% setup wdir indices
% could fill NaNs with EWS data, but actual wdir is different from PWS ~70%
% of the time (even though the wind rose looks similar)
ind_N=bee_dataset.N_SE_rest==1;
ind_SE=bee_dataset.N_SE_rest==2;
ind_rest=bee_dataset.N_SE_rest==3;

% load data grouped for flexpart
load('/home/kristof/work/documents/paper_bro/data/BEE_dataset_flexpart.mat')

% 2015 measurements exluded later
    
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
    
    eval(['c_list=' c_scale '(4);'])
    
    plot_data=bee_dataset;
    
    bad_ind=find(bee_dataset.times.Year==2017 & bee_dataset.times.Month==4);
    plot_data.bro_col(bad_ind)=NaN;
    plot_data.aer_ext(bad_ind)=NaN;
    plot_data.aer_halfmicron(bad_ind)=NaN;
    plot_data.o3_surf(bad_ind)=NaN;
    
    plot_data=[plot_data;plot_data(end,:)];
    plot_data.times(end)=datetime(2017,05,13);
    plot_data.bro_col(end)=NaN;
    plot_data.aer_ext(end)=NaN;
    plot_data.aer_halfmicron(end)=NaN;
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
               
    ylabel('BrO VCD (molec\timescm^{-2})')
    
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

    ylabel('D_p > 0.5 \mum (cm^{-1})')
    
    text(labelx,labely,'d)','color','k','FontWeight','bold','Units','normalized')    
    
    %%% plot formatting
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
%     ll.FontSize=11;

    save_pdf(save_figs, '2x2_box')
    
    
    
end

%% daily mean BrO, with min-max range
if bro_dailymean
    
    % check if daily mean BrO was saved already
    if ~exist('/home/kristof/work/documents/paper_bro/data/daily_mean_BrO.mat','file')

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
        
        save('/home/kristof/work/documents/paper_bro/data/daily_mean_BrO.mat','dmean_bro')
        
    else
        % if already saved 
        load('/home/kristof/work/documents/paper_bro/data/daily_mean_BrO.mat')
    end
    
    % plot figure
    figure
    set(gcf, 'Position', [100, 100, 1000, 600]);
    fig_ax = tight_subplot(4,1,[0.065,0.04],[0.12,0.068],[0.104,0.04]);
    
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
                      'facecolor','b','edgecolor','b')
            rectangle('position',...
                      [plot_x(dy)-box_w/2,plot_mean(dy),box_w,1e12],...
                      'facecolor','k','edgecolor','k')
                  
        end

        grid on
        ylim([-5e12,16e13])
        xlim([datenum(yr,3,1),datenum(yr,6,2)])
        
        set(gca, 'XTick', datenum([[yr,3,1];[yr,3,15];[yr,4,1];[yr,4,15];...
                                   [yr,5,1];[yr,5,15];[yr,6,1]]))
        set(gca, 'XTicklabel', {'01/03','15/03','01/04','15/04','01/05','15/05','01/06'})
        
        if yr~=2019
            set(gca,'xticklabel',[])
        else
            xlabel('Date (UTC)')
        end
        
        if yr==2017
            ylb=ylabel('BrO VCD (molec/cm^2)'); 
            ylb.Position(1)=ylb.Position(1)-3;
            ylb.Position(2)=ylb.Position(2)-max(dmean_bro.mean);
        end
        
        text(0.9,0.8,num2str(yr), 'color','k','Units','normalized',...
        'fontsize',fig_fs,'fontweight','bold')
        
    end
    
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)

    save_pdf(save_figs, 'BrO_dmean')
    
end


%% plot relationship to wind speed
if weather_corr
    
    edges=0:18;
    
    txt_pos=0.92;
         
    
    % set up figure
    figure
    set(gcf, 'Position', [100, 100, 1000, 600]);
    fig_ax = tight_subplot(2,3,[0.07,0.07],[0.08,0.1],[0.08,0.03]);

    %%%%%% wind speed
    % northerly winds, mean: 354 deg from gaussian fit, +-30 deg
    axes(fig_ax(1))

    dscatter(bee_dataset.wspd_ms(ind_N), bee_dataset.bro_col(ind_N))
    hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind_N),bee_dataset.bro_col(ind_N),edges)
    
    text(0.95,txt_pos,'{\bfa)} Wind: N', 'color','k','Units','normalized',...
        'fontsize',fig_fs,'HorizontalAlignment','right')
    ylabel('BrO VCD (molec/cm^2)')
    xlabel('Wind speed (m/s)')
    ylim([0,12]*1e13)
    xlim([0,15])   
    
    % southeasterly winds, mean: 123 deg from gaussian fit, +-30 deg  
    axes(fig_ax(2))

    dscatter(bee_dataset.wspd_ms(ind_SE), bee_dataset.bro_col(ind_SE))
    hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind_SE),bee_dataset.bro_col(ind_SE),edges)
    
    text(0.95,txt_pos,'{\bfb)} Wind: SE', 'color','k','Units','normalized',...
        'fontsize',fig_fs,'HorizontalAlignment','right')
    ylabel('BrO VCD (molec/cm^2)')
    xlabel('Wind speed (m/s)')
    ylim([0,12]*1e13)
    xlim([0,15])   

    % everything else
    axes(fig_ax(3))

    dscatter(bee_dataset.wspd_ms(ind_rest), bee_dataset.bro_col(ind_rest))
    hold on, box on
    plot_mean_std(bee_dataset.wspd_ms(ind_rest),bee_dataset.bro_col(ind_rest),edges)
    
    text(0.95,txt_pos,'{\bfc)} Wind: other', 'color','k','Units','normalized',...
        'fontsize',fig_fs,'HorizontalAlignment','right')
    ylabel('BrO VCD (molec/cm^2)')
    xlabel('Wind speed (m/s)')
    ylim([0,12]*1e13)
    xlim([0,15])   
    


    
    
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)

    save_pdf(save_figs, 'BrO_weather')
    
end


%% back trajectories and surface sensitivity

if sens_meanBrO
    
    % load sensitivity and back trajectory data
    load('/home/kristof/work/documents/paper_bro/data/BrO_back_runs_v1_3day.mat')
        
    figure
    set(gcf, 'Position', [100, 100, 1000, 800]);
    fig_ax = tight_subplot(2,3,[0.08,0.02],[0.04,0.1],[0.04,0.04]);
    
    % N winds and above average BrO cols
    axes(fig_ax(1))
    plot_ind=(bee_fp.wdir==1 & bee_fp.bro_m==2 & bee_fp.mean_time.Year>2015); 
    
    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind);
    title('{\bfa)} Wind: N','fontweight','normal','FontSize',fig_fs)
    
    % SE winds and above average BrO cols
    axes(fig_ax(2))
    plot_ind=(bee_fp.wdir==2 & bee_fp.bro_m==2 & bee_fp.mean_time.Year>2015); 
    
    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind);
    title('{\bfb)} Wind: SE','fontweight','normal','FontSize',fig_fs)

    % other winds and above average BrO cols
    axes(fig_ax(3))
    plot_ind=(bee_fp.wdir>=3 & bee_fp.bro_m==2 & bee_fp.mean_time.Year>2015);
    
    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind);
    title('{\bfc)} Wind: other','fontweight','normal','FontSize',fig_fs)
    
    %%%%%%%%%%%%%%%%%
    
    % N winds and ?
    axes(fig_ax(4))
    plot_ind=(bee_fp.wdir==1 & bee_fp.bro_pc==10 & bee_fp.mean_time.Year>2015); 
    
    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind);
    title('{\bfd)} Wind: N','fontweight','normal','FontSize',fig_fs)
    
    % SE winds and ?
    axes(fig_ax(5))
    plot_ind=(bee_fp.wdir==2 & bee_fp.bro_pc==10 & bee_fp.mean_time.Year>2015); 
    
    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind);
    title('{\bfe)} Wind: SE','fontweight','normal','FontSize',fig_fs)

    % other winds ?
    axes(fig_ax(6))
    plot_ind=(bee_fp.wdir>=3 & bee_fp.bro_pc==10 & bee_fp.mean_time.Year>2015);
    
    plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind);
    title('{\bff)} Wind: other','fontweight','normal','FontSize',fig_fs)
    
    if save_figs, save_pdf(2, 'sens_map'), end
    
end



%% correlation with coarse mode particles
if plot_ssa
    
    
    figure
    set(gcf, 'Position', [100, 100, 1000, 700]);
    fig_ax = tight_subplot(2,2,[0.08,0.05],[0.12,0.07],[0.11,0.05]);
    
    txt_pos=0.92;
    ind_ssa=(bee_dataset.aer_halfmicron <= 100);

    plot_type='hm';
    
    switch plot_type
        case 'sm'
            plot_data=bee_dataset.aer_supermicron;
            xlim_end=1;
            x_label='D_P > 1 \mum (cm^{-3})';
        case 'hm'
            plot_data=bee_dataset.aer_halfmicron;
            xlim_end=6;
            x_label='D_P > 0.5 \mum (cm^{-3})';
    end
    
    % all aer data
    axes(fig_ax(1))
    dscatter(plot_data(ind_ssa), bee_dataset.bro_col(ind_ssa)), hold on, box on
    
    plot_fit_line(plot_data(ind_ssa),bee_dataset.bro_col(ind_ssa),fig_fs)
    
    ylim([0,16]*1e13)
    xlim([0,xlim_end])    
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'{\bfa)} Wind: all', 'color','k','Units','normalized',...
        'fontsize',fig_fs,'HorizontalAlignment','right')
    
    % Northerly winds only
    axes(fig_ax(2))
    ind=(ind_ssa & ind_N);
    dscatter(plot_data(ind), bee_dataset.bro_col(ind)), hold on, box on
    
    plot_fit_line(plot_data(ind),bee_dataset.bro_col(ind),fig_fs)    
    
    ylim([0,16]*1e13)
    xlim([0,xlim_end])    

    text(0.95,txt_pos,'{\bfb)} Wind: N', 'color','k','Units','normalized',...
        'fontsize',fig_fs,'HorizontalAlignment','right')
    
    % Southeasterly winds only
    axes(fig_ax(3))
    ind=(ind_ssa & ind_SE);
    dscatter(plot_data(ind), bee_dataset.bro_col(ind)), hold on, box on
    
    plot_fit_line(plot_data(ind),bee_dataset.bro_col(ind),fig_fs)    
    
    ylim([0,16]*1e13)
    xlim([0,xlim_end])    
    xlabel(x_label)
    ylabel('BrO VCD (molec/cm^2)')

    text(0.95,txt_pos,'{\bfc)} Wind: SE', 'color','k','Units','normalized',...
        'fontsize',fig_fs,'HorizontalAlignment','right')
    
    % all other wind directions
    axes(fig_ax(4))
    ind=(ind_ssa & ind_rest);
    dscatter(plot_data(ind), bee_dataset.bro_col(ind)), hold on, box on
    
    plot_fit_line(plot_data(ind),bee_dataset.bro_col(ind),fig_fs)    
    
    ylim([0,16]*1e13)
    xlim([0,xlim_end])    
    xlabel(x_label)

    text(0.95,txt_pos,'{\bfd)} Wind: other', 'color','k','Units','normalized',...
        'fontsize',fig_fs,'HorizontalAlignment','right')
    
    
    set(findall(gcf,'-property','FontSize'),'FontSize',fig_fs)
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
    save_pdf(save_figs, 'ssa_corr')
    
end

%% SI contact as a function of BrO percentile
% using mean BrO column for each trajectory
if si_contact
    
    eval(['c_list=' c_scale '(10);'])
%     c_list_tmp=flipud(c_list_tmp);

    figure
    set(gcf, 'Position', [100, 100, 1000, 600]);
    fig_ax = tight_subplot(2,1,[0.06,0.01],[0.085,0.06],[0.086,0.17]);
    % fig_ax = tight_subplot(2,1,[0.07,0.01],[0.08,0.06],[0.08,0.16]);

    % select data to plot
    plot_type_all={'FYSI','MYSI'};

    for sp=1:2

        axes(fig_ax(sp))

        plot_type=plot_type_all{sp};

        switch plot_type
            case 'FYSI'
                box_data=bee_fp.fysi_3day;
            case 'MYSI'
                box_data=bee_fp.mysi_3day;
            case 'water'
                box_data=bee_fp.water_3day;
            case 'land'
                box_data=bee_fp.land_3day;
        end

        % create indices to plot bar plots with
        % each number in ascending order plots a different bar, code below sorts
        % BrO percentiles by wind direction, since wdir indices are 1-4  
        % (gives 100+ for N, 200+ for SE, 300+ for other)
        bee_fp.wdir(bee_fp.wdir==4)=3; % merge undetermined wdirs with 'other'
        box_group=bee_fp.bro_pc+100*bee_fp.wdir;
        % box_group=bee_fp.ssa_pc+100*bee_fp.wdir;

        % set positions to separate bar groups by wdir
        a=0.08; % controls width of each group
        nbars=10; % number of bars in each group (number of unique bro_pc indices)
        box_pos=sort([1:a:1+a*(nbars-1),2:a:2+a*(nbars-1),3:a:3+a*(nbars-1)]);
        
        % make box plot
        boxplot(box_data,box_group,'positions',box_pos,'colors',c_list,...
            'jitter',0.5,'symbol','.')

        hold on, box on
        
        % add mean as well
        group_ind=unique(box_group);
        c_ind=repmat(1:nbars,1,length(group_ind)/nbars);
        for i=1:length(group_ind)
            plot(box_pos(i), nanmean(box_data(box_group==group_ind(i))),'x','color',c_list(c_ind(i),:))
        end
        
        % set xtick position to mean of each group
        tick_pos=[];
        for i=1:3
            tick_pos=[tick_pos,mean(box_pos((nbars)*(i-1)+1:(nbars)*i))];
        end
        
        set(findobj(gca,'type','line'),'linew',box_lw)
        set(findall(gca,'tag','Outliers'),'MarkerSize',box_outlier);
        
        set(gca,'xtick',tick_pos)
        set(gca,'xticklabel',{'N','SE','other'})

        xlim([box_pos(1)-0.15,box_pos(end)+0.15])        
        ylim([-0.7,10.7]*1e13)
        set(gca,'ygrid','on')

        set(findall(gca,'-property','FontSize'),'FontSize',fig_fs)
        
        if sp==1

            ll=legend(flipud(findall(gca,'Tag','Box')),...
                   {'<10^t^h','10^t^h-20^t^h','20^t^h-30^t^h','30^t^h-40^t^h','40^t^h-50^t^h',...
                    '50^t^h-60^t^h','60^t^h-70^t^h','70^t^h-80^t^h','80^t^h-90^t^h','>90^t^h'},...
                   'Orientation','vertical','location','east');

            ll.FontSize=11;

            % [left bottom width height]
    %         ll_pos=get(ll,'position');
            ll.Position(1)=0.845;
            ll.Position(2)=0.3;

            text(1.095,0.48,...
                 sprintf('BrO column\npercentiles'), 'color','k','Units','normalized',...
                 'fontsize',13,'HorizontalAlignment','center')

        else
            xlb=xlabel('Mean wind direction','fontsize',fig_fs);
%             xlb.Position(2)=xlb.Position(2)*1.9;
        end

        ylabel([plot_type ' contact (s m^2)'])

    end
    
    set(findall(gcf,'-property','FontName'),'FontName',fig_font)
    
    save_pdf(save_figs, 'SI_contact')
    
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
            print(h,f_out,'-djpeg','-r400','-opengl')
%             % png images
%             f_out=['/home/kristof/work/documents/paper_sat-val/figures/' fname '.png'];
%             print(h,f_out,'-dpng','-r400','-opengl')
            
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

function plot_mean_std(xx,yy,edges)

    tmp_mean=[];
    tmp_std=[];
    
    for i=1:length(edges)-1
        tmp_mean=[tmp_mean,nanmean(yy(xx>=edges(i) & xx < edges(i+1)))];
        tmp_std=[tmp_std,nanstd(yy(xx>=edges(i) & xx < edges(i+1)))];
    end    
    
    plot_x=edges(1:end-1)+diff(edges)/2;
    
    plot(plot_x,tmp_mean,'k-','linewidth',2)
    plot(plot_x,tmp_mean+tmp_std,'k--','linewidth',2)
    plot(plot_x,tmp_mean-tmp_std,'k--','linewidth',2)   

end

function plot_fit_line(xx,yy,text_size,loc)

    if nargin==3, loc='left'; end

    [slope, y_int, R2] = line_fit(xx,yy);
    plot([min(xx):max(xx)/20:max(xx)+2],...
          slope(1)*[min(xx):max(xx)/20:max(xx)+2]+y_int(1),'linewidth',1)

    if strcmp(loc,'left')
        text(0.05,0.92,sprintf('N=%.0f', size(xx,1)), 'color','k','Units','normalized',...
            'fontsize',text_size)

        text(0.05,0.76,sprintf('{R^2}=%.2f', R2),'color','k','Units','normalized',...
            'fontsize',text_size)

%         text(0.05,0.85,sprintf('m=%.2g', slope(1)),'color','k','Units','normalized',...
%             'fontsize',text_size)

    elseif strcmp(loc,'right')
        text(0.95,0.95,sprintf('N=%.0f', size(xx,1)), 'color','k','Units','normalized',...
            'fontsize',text_size,'HorizontalAlignment','right')

        text(0.95,0.90,sprintf('{R^2}=%.2f', R2),'color','k','Units','normalized',...
            'fontsize',text_size,'HorizontalAlignment','right')

        text(0.95,0.85,sprintf('m=%.2g', slope(1)),'color','k','Units','normalized',...
            'fontsize',text_size,'HorizontalAlignment','right')
    end
    
end

function plot_sens_map(latitude,longitude,sensitivities,trajectories,plot_ind)

    if length(plot_ind)~=size(sensitivities,3)
        error('plot_ind must be a boolean array, same length as number of flexpart runs');
    end

    load coast;
    ax = worldmap([61,90], [-180,180]);
    geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])
    hold on

    sensitivities(sensitivities==0)=NaN;

    sens_tmp=nanmean(sensitivities(:,:,plot_ind),3);

    sens_tmp(sens_tmp<0.1)=0.1;

    % plot on log scale
    sens_tmp=log10(sens_tmp);
    % replace max log value with next largest integer if difference is <0.3
    [max_val,max_ind] = max(sens_tmp(:));
    if max_val-floor(max_val) >0.7
        sens_tmp(max_ind)=ceil(max_val);
        max_val=ceil(max_val);
    end

    % plot sensitivity
    surfm(latitude,longitude,sens_tmp,'facecolor', 'interp','facealpha',0.75);
    colormap(flipud(hot))

    %     if plot_subplot && i==num_plots
    %         % add colorbar with manual labels (convert back to lin space)
    %         cb=colorbar('position',[0.873,0.29,0.026,0.47]);
    %         cb_lim=-1:ceil(max_val);
    %         set(cb,'YTick',cb_lim)
    %         set(cb,'YTick',cb_lim,'YTickLabel',cellstr(num2str(power(10,cb_lim)')))
    %         ylabel(cb,'Mean surf. sens. (s)','FontSize',12)
    %     end

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
