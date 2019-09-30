function bro_corr_plots
% generate correlation plots of BrO columns with other parameters

% select plots
wspd_corr=0;
wdir_corr=0;
t_rl_corr=0;
opc_corr=0;
smps_corr=0;
o3_corr=0;

ptom_comp=1;

plot_availability=0;

text_size=11;
    
% true for using 0-4 km columns in figures
% set to false to use ratio of part col below lab to total column (set part col alt below)
plot_column=1;

% load data
load('/home/kristof/work/BEEs/BEE_dataset_all.mat');

% setup plotting indices
ind_N=bee_dataset.N_SE_rest==1;
ind_SE=bee_dataset.N_SE_rest==2;
ind_rest=bee_dataset.N_SE_rest==3;

ind_t_rl=~isnan(bee_dataset.T_PWS);
ind_opc=(bee_dataset.OPC_supermicron <= 1);
ind_smps=~isnan(bee_dataset.SMPS_100_500);
ind_o3=~isnan(bee_dataset.o3_surf);

%% set up plot variable
if plot_column
    plot_var=bee_dataset.bro_col;
else % plot ratio
    plot_var=bee_dataset.bro_col_ratio;   
end


%% Plot correlations

if wspd_corr
    
% % %     WindRose(data.WindDir,data.WindSpd,'AngleNorth',0,'AngleEast',90,...
% % %              'FreqLabelAngle',45,'vWinds',[0 3 6 9 12 15]);
    
    if plot_column, txt_pos=0.92; else txt_pos=0.08; end
         
    % all winds
    figure
    subplot(221), hold on, box on
    ind=~isnan(bee_dataset.wspd_ms);
    dscatter(bee_dataset.wspd_ms(ind), plot_var(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),plot_var(ind),plot_column)

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    if plot_column, ylim([0,12]*1e13), else ylim([0.28,1]), end
    xlim([0,15])   
    if plot_column
        ylabel('BrO VCD_{0-4 km} (molec/cm^2)')
    else
         ylabel('VCD_{0-0.6 km} / VCD_{0-4 km}')
    end
    
    % northerly winds, mean: 354 deg from gaussian fit, +-30 deg
    subplot(222), hold on, box on

    dscatter(bee_dataset.wspd_ms(ind_N), plot_var(ind_N))
    plot_mean_std(bee_dataset.wspd_ms(ind_N),plot_var(ind_N),plot_column)
    
    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    if plot_column, ylim([0,12]*1e13), else ylim([0.28,1]), end
    xlim([0,15])   
    
    % southeasterly winds, mean: 123 deg from gaussian fit, +-30 deg  
    subplot(223), hold on, box on

    dscatter(bee_dataset.wspd_ms(ind_SE), plot_var(ind_SE))
    plot_mean_std(bee_dataset.wspd_ms(ind_SE),plot_var(ind_SE),plot_column)
    
    text(0.95,txt_pos,'Wind: 123° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    if plot_column, ylim([0,12]*1e13), else ylim([0.28,1]), end
    xlim([0,15])    
    xlabel('Wind speed (m/s)')
    if plot_column
        ylabel('BrO VCD_{0-4 km} (molec/cm^2)')
    else
         ylabel('VCD_{0-0.6 km} / VCD_{0-4 km}')
    end

    % everything else
    subplot(224), hold on, box on

    dscatter(bee_dataset.wspd_ms(ind_rest), plot_var(ind_rest))
    plot_mean_std(bee_dataset.wspd_ms(ind_rest),plot_var(ind_rest),plot_column)
    
    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    if plot_column, ylim([0,12]*1e13), else ylim([0.28,1]), end
    xlim([0,15])    
    xlabel('Wind speed (m/s)')
    
    
end

if wdir_corr
    
    figure
    ind=find(~isnan(bee_dataset.wdir) & bee_dataset.wdir>0);
    dscatter(bee_dataset.wdir(ind), plot_var(ind))
    
end

if t_rl_corr
    
%     figure
%     dscatter(bee_dataset.T_PWS(ind), plot_var(ind))
    
    figure
    txt_pos=0.92;

    % all data
    subplot(221), hold on, box on
    dscatter(bee_dataset.T_PWS(ind_t_rl), plot_var(ind_t_rl))
    
    plot_fit_line(bee_dataset.T_PWS(ind_t_rl),plot_var(ind_t_rl),text_size)
    
    ylim([0,16]*1e13)
    xlim([-45,5])    
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Northerly winds only
    txt_pos=0.92;
    subplot(222), hold on, box on
    ind=(ind_t_rl & ind_N);
    dscatter(bee_dataset.T_PWS(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.T_PWS(ind),plot_var(ind),text_size)    
    
    ylim([0,16]*1e13)
    xlim([-45,5])    

    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Southeasterly winds only
    subplot(223), hold on, box on
    ind=(ind_t_rl & ind_SE);
    dscatter(bee_dataset.T_PWS(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.T_PWS(ind),plot_var(ind),text_size)    
    
    ylim([0,16]*1e13)
    xlim([-45,5])    
    xlabel('PWS T (°C)')
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'Wind: 123° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % all other wind directions
    subplot(224), hold on, box on
    ind=(ind_t_rl & ind_rest);
    dscatter(bee_dataset.T_PWS(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.T_PWS(ind),plot_var(ind),text_size)    
    
    ylim([0,16]*1e13)
    xlim([-45,5])    
    xlabel('PWS T (°C)')

    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    
end

if opc_corr

    figure
    txt_pos=0.08;

    % all OPC data
    subplot(221), hold on, box on
    dscatter(bee_dataset.OPC_supermicron(ind_opc), plot_var(ind_opc))
    
    plot_fit_line(bee_dataset.OPC_supermicron(ind_opc),plot_var(ind_opc),text_size)
    
    ylim([0,16]*1e13)
    xlim([0,1])    
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Northerly winds only
    txt_pos=0.92;
    subplot(222), hold on, box on
    ind=(ind_opc & ind_N);
    dscatter(bee_dataset.OPC_supermicron(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.OPC_supermicron(ind),plot_var(ind),text_size)    
    
    ylim([0,16]*1e13)
    xlim([0,1])    

    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Southeasterly winds only
    subplot(223), hold on, box on
    ind=(ind_opc & ind_SE);
    dscatter(bee_dataset.OPC_supermicron(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.OPC_supermicron(ind),plot_var(ind),text_size)    
    
    ylim([0,16]*1e13)
    xlim([0,1])    
    xlabel('OPC 1-10 \mum (cm^{-3})')
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'Wind: 123° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % all other wind directions
    subplot(224), hold on, box on
    ind=(ind_opc & ind_rest);
    dscatter(bee_dataset.OPC_supermicron(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.OPC_supermicron(ind),plot_var(ind),text_size)    
    
    ylim([0,16]*1e13)
    xlim([0,1])    
    xlabel('OPC 1-10 \mum (cm^{-3})')

    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
end

if smps_corr
    
    figure, hold on

    subplot(221), hold on, box on
    dscatter(bee_dataset.SMPS_100_500(ind_smps), plot_var(ind_smps))
    plot_fit_line(bee_dataset.SMPS_100_500(ind_smps),plot_var(ind_smps),text_size)
    xlim([0,500])
    ylim([0,16]*1e13)
    
    subplot(222), hold on, box on
    ind=(ind_smps & ind_N);
    dscatter(bee_dataset.SMPS_100_500(ind), plot_var(ind))
    plot_fit_line(bee_dataset.SMPS_100_500(ind),plot_var(ind),text_size)
    xlim([0,500])
    ylim([0,13]*1e13)

    subplot(223), hold on, box on
    ind=(ind_smps & ind_SE);
    dscatter(bee_dataset.SMPS_100_500(ind), plot_var(ind))
    plot_fit_line(bee_dataset.SMPS_100_500(ind),plot_var(ind),text_size)
    xlim([0,500])
    ylim([0,13]*1e13)

    subplot(224), hold on, box on
    ind=(ind_smps & ind_rest);
    dscatter(bee_dataset.SMPS_100_500(ind), plot_var(ind))
    plot_fit_line(bee_dataset.SMPS_100_500(ind),plot_var(ind),text_size)
    xlim([0,500])
    ylim([0,13]*1e13)
    
end



if o3_corr
    
%     figure, hold on
%     dscatter(bee_dataset.surf_o3(ind_o3), plot_var(ind_o3))
% %     plot_fit_line(bee_dataset.surf_o3(ind_o3),plot_var(ind_o3),text_size,'right')
%     ylim([0,9]*1e13)
%     
%     figure, hold on
%     ind=(ind_o3, & ind_N);
%     dscatter(bee_dataset.surf_o3(ind), plot_var(ind))
% %     plot_fit_line(bee_dataset.surf_o3(ind),plot_var(ind),text_size,'right')
%     ylim([0,9]*1e13)
% 
%     figure, hold on
%     ind=(ind_o3 & ind_SE);
%     dscatter(bee_dataset.surf_o3(ind), plot_var(ind))
% %     plot_fit_line(bee_dataset.surf_o3(ind),plot_var(ind),text_size,'right')
%     ylim([0,9]*1e13)
% 
%     figure, hold on
%     ind=(ind_o3 & ind_rest);
%     dscatter(bee_dataset.surf_o3(ind), plot_var(ind))
% %     plot_fit_line(bee_dataset.surf_o3(ind),plot_var(ind),text_size,'right')
%     ylim([0,9]*1e13)
    

    % bar plot of BrO in ozone bins (ignores NaN)
    
    bro_limit=mean(bee_dataset.bro_col);
%     bro_sort=sort(bee_dataset.bro_col);
%     bro_limit=bro_sort(floor(length(bee_dataset.bro_col)*0.75)+1);
    
    o3_0_10=(bee_dataset.surf_o3<10 & bee_dataset.bro_col>=bro_limit);
    o3_10_20=(bee_dataset.surf_o3>=10 & bee_dataset.surf_o3<20 & ...
              bee_dataset.bro_col>=bro_limit);
    o3_20_30=(bee_dataset.surf_o3>=20 & bee_dataset.surf_o3<30 & ...
              bee_dataset.bro_col>=bro_limit);
    o3_30_up=(bee_dataset.surf_o3>=30 & bee_dataset.bro_col>=bro_limit);
    
    by_wind=[ [sum(o3_0_10 & ind_N),...
                sum(o3_0_10 & ind_SE),...
                sum(o3_0_10 & ind_rest) ];...
              [sum(o3_10_20 & ind_N),...
                sum(o3_10_20 & ind_SE),...
                sum(o3_10_20 & ind_rest) ];...
              [sum(o3_20_30 & ind_N),...
                sum(o3_20_30 & ind_SE),...
                sum(o3_20_30 & ind_rest) ];...
              [sum(o3_30_up & ind_N),...
                sum(o3_30_up & ind_SE),...
                sum(o3_30_up & ind_rest) ] ];
            
    figure
    bar(by_wind./sum(sum(by_wind)),'stacked');
    legend('354° \pm 30°','123° \pm 30°','Other','location','north')
    xlabel('Surface O_3 (ppbv)')
    ylabel('p for BrO VCD_{0-4 km} > mean')
    xlim([0.5,4.5])
    
    c = {'0-10','10-20','20-30','>30'};
    set(gca,'xticklabel',c)
    
    
end


if ptom_comp
    
    ms=8;
    figure
    
    load('/home/kristof/work/BEEs/ptom_pcr.mat');
    
    % 2015 data
    subplot(221), hold on, box on, grid on
    ind=year(bee_dataset.times)==2015;
    plot(bee_dataset.times(ind), bee_dataset.bro_col(ind),'ks','markersize',ms,'linewidth',1)
    plot(bee_dataset.times(ind), bee_dataset.bro_col_ptom(ind),'ro','markersize',ms,'linewidth',1)
    plot(bee_dataset.times(ind), pca_bro(ind),'b.','markersize',ms+2,'linewidth',1)
    legend('MAX-DOAS','pTOMCAT','PCR with pTOMCAT')
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')
    ylim([0,3]*1e14)
    xlim([datenum('2015/02/27'),datenum('2015/06/03')])
    set(gca,'XTick',[datenum('2015/03/01'),datenum('2015/04/01'),...
                     datenum('2015/05/01'),datenum('2015/06/01')])   
    
    text(0.04,0.90,'2015', 'color','k','Units','normalized',...
        'fontsize',text_size,'fontweight','bold')
    
    % 2016 data
    subplot(222), hold on, box on, grid on
    ind=year(bee_dataset.times)==2016;
    plot(bee_dataset.times(ind), bee_dataset.bro_col(ind),'ks','markersize',ms,'linewidth',1)
    plot(bee_dataset.times(ind), bee_dataset.bro_col_ptom(ind),'ro','markersize',ms,'linewidth',1)
    plot(bee_dataset.times(ind), pca_bro(ind),'b.','markersize',ms+2,'linewidth',1)
    ylim([0,3]*1e14)
    xlim([datenum('2016/02/28 07:00'),datenum('2016/05/13')])
    set(gca,'XTick',[datenum('2016/03/01'),datenum('2016/04/01'),datenum('2016/05/01')])
    
    text(0.96,0.90,'OPC available', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    text(0.04,0.90,'2016', 'color','k','Units','normalized',...
        'fontsize',text_size,'fontweight','bold')
    
    % 2017 data
    subplot(223), hold on, box on, grid on
    ind=year(bee_dataset.times)==2017;
    plot(bee_dataset.times(ind), bee_dataset.bro_col(ind),'ks','markersize',ms,'linewidth',1)
    plot(bee_dataset.times(ind), bee_dataset.bro_col_ptom(ind),'ro','markersize',ms,'linewidth',1)
    plot(bee_dataset.times(ind), pca_bro(ind),'b.','markersize',ms+2,'linewidth',1)
    ylim([0,3]*1e14)
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')
    xlim([datenum('2017/02/28 5:00'),datenum('2017/04/06')])
    set(gca,'XTick',[datenum('2017/03/01'),datenum('2017/04/01')])

    text(0.96,0.88,'OPC, O_3 available', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    text(0.04,0.90,'2017', 'color','k','Units','normalized',...
        'fontsize',text_size,'fontweight','bold')
    
    % 2018 data
    subplot(224), hold on, box on, grid on
    ind=year(bee_dataset.times)==2018;
    plot(bee_dataset.times(ind), bee_dataset.bro_col(ind),'ks','markersize',ms,'linewidth',1)
    plot(bee_dataset.times(ind), bee_dataset.bro_col_ptom(ind),'ro','markersize',ms,'linewidth',1)
    plot(bee_dataset.times(ind), pca_bro(ind),'b.','markersize',ms+2,'linewidth',1)
    ylim([0,3]*1e14)
    xlim([datenum('2018/02/27'),datenum('2018/06/03')])
    set(gca,'XTick',[datenum('2018/03/01'),datenum('2018/04/01'),...
                     datenum('2018/05/01'),datenum('2018/06/01')])   

    text(0.96,0.88,'OPC, O_3 available', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    text(0.04,0.90,'2018', 'color','k','Units','normalized',...
        'fontsize',text_size,'fontweight','bold')
    
    
%     figure, hold on
%     plot([0,2.55e14],[0,2.55*1e14],'k--')
%     dscatter(bee_dataset.bro_col(coincidences(:,1)),ptom_col(coincidences(:,2)))
%     plot_fit_line(bee_dataset.bro_col(coincidences(:,1)),ptom_col(coincidences(:,2)),text_size,'right')
%     xlim([0,2.55]*1e14)
%     ylim([0,2.55]*1e14)
%     pbaspect([1 1 1])
    
end

%% plot data availability
% consider March-May only
if plot_availability
    
    figure
    fig_ax = tight_subplot(2,2,[0.1,0.07],[0.15,0.11],[0.07,0.04]);
    
    % BrO
    plot_avail(times, 1, fig_ax)   
    title('BrO')
    
    % OPC
    plot_avail(opc_time, 2, fig_ax)
    title('OPC')
    
    % SMPS
%     tmp=smps_time;
%     tmp=dateshift(tmp,'start','hour','nearest');
%     tmp=unique(tmp);
%     
%     plot_avail(tmp, 3)        

    % Surf_ozone
    tmp=surf_o3.DateTime;
    tmp=dateshift(tmp,'start','hour','nearest');
    tmp=unique(tmp);
    
    plot_avail(tmp, 3, fig_ax)   
    title('Surface ozone')

    % PWS
    tmp=data.DateTime;
    tmp=dateshift(tmp,'start','hour','nearest');
    tmp=unique(tmp);
    
    plot_avail(tmp, 4, fig_ax)
    title('PEARL weather station')

    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(findall(gcf,'-property','FontName'),'FontName','Arial') 
    set(gcf, 'Position', [100, 100, 1000, 300]);

end
end

%% extra functions
function plot_avail(tmp, plotnum, fig_ax)

    tmp(month(tmp)>5 | month(tmp)<3)=[];
    yr=year(tmp);
    tmp.Year=0;

    xt=[61, 75, 92, 106, 122, 136, 152];
    yt=2015:2018;
    
%     subplot(2,2,plotnum)
    axes(fig_ax(plotnum))
    plot(tmp,yr, 'marker','+','linestyle','none')
    set(gca,'XTick',xt)
    set(gca,'YTick',yt)
    set(gca,'yticklabel',{'2015','','','2018'})
    set(gca,'XGrid','on')
    set(gca,'GridAlpha',.5)
    ylim([2014.5,2018.5])
    xlim([60,154])

    if any(plotnum==[1,2]), set(gca,'xticklabel',{[]}); end


end

function plot_mean_std(xx,yy,plot_column)

    if nargin==2, plot_column=1; end
    
    tmp_mean=[];
    tmp_std=[];
    
    for i=0:17
        tmp_mean=[tmp_mean,nanmean(yy(xx>=i & xx < i+1))];
        tmp_std=[tmp_std,nanstd(yy(xx>=i & xx < i+1))];
    end    
    
    
    if ~plot_column, plot([0,17],[0.83,0.83],'b-','linewidth',1), end

    plot([0:17]+0.5,tmp_mean,'k-','linewidth',2)
    plot([0:17]+0.5,tmp_mean+tmp_std,'k--','linewidth',2)
    plot([0:17]+0.5,tmp_mean-tmp_std,'k--','linewidth',2)   

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





