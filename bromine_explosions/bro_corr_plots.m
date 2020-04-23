function [r_all,lag_all]=bro_corr_plots()
% generate correlation plots of BrO columns with other parameters

% select plots
windrose=0;
wspd_corr=0;
wdir_corr=0;
T_rl_corr=0;
T_ews_corr=0;
sonde_dT_corr=0;
sonde_dT_o3=0;
ssa_corr=0;
ssa_corr_oneplot=1;
smps_corr=0;
o3_corr=0;

aod_corr=0;
aod_to_aer=0;

aod_wspd=0;
aer_wspd=0;
smps_wspd=0;

SI_corr=0;
SI_corr_flip=0;

SI_to_SI=0;

ptom_comp=0;

plot_availability=0;

text_size=11;
    
% true for using 0-4 km columns in figures
% set to false to use ratio of part col below lab to total column (set part col alt below)
plot_column=1;

% load data
load('/home/kristof/work/BEEs/BEE_dataset_all.mat');
bee_dataset(bee_dataset.times.Year==2015,:)=[];
% bee_dataset(bee_dataset.bro_col<prctile(bee_dataset.bro_col,90),:)=[];
% bee_dataset(bee_dataset.sonde_dT<7,:)=[];


% setup plotting indices
ind_N=bee_dataset.N_SE_rest==1;
ind_SE=bee_dataset.N_SE_rest==2;
ind_rest=bee_dataset.N_SE_rest==3;

ind_t_rl=~isnan(bee_dataset.T_PWS);
ind_t_ews=~isnan(bee_dataset.T_EWS);
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
    
    edges=0:18;
    
    if plot_column, txt_pos=0.92; else txt_pos=0.08; end
         
    % all winds
    figure
    subplot(221), hold on, box on
    ind=(~isnan(bee_dataset.wspd_ms_EWS) & ~isnan(plot_var));
    dscatter(bee_dataset.wspd_ms_EWS(ind), plot_var(ind))
    plot_mean_std(bee_dataset.wspd_ms_EWS(ind),plot_var(ind),edges,plot_column)

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
%     if plot_column, ylim([0,12]*1e13), else ylim([0.28,1]), end
    xlim([0,15])   
    if plot_column
        ylabel('BrO VCD_{0-4 km} (molec/cm^2)')
    else
         ylabel('VCD_{0-0.6 km} / VCD_{0-4 km}')
    end
    
    % northerly winds, mean: 354 deg from gaussian fit, +-30 deg
    subplot(222), hold on, box on

    dscatter(bee_dataset.wspd_ms_EWS(ind_N), plot_var(ind_N))
    plot_mean_std(bee_dataset.wspd_ms_EWS(ind_N),plot_var(ind_N),edges,plot_column)
    
    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    if plot_column, ylim([0,12]*1e13), else ylim([0.28,1]), end
    xlim([0,15])   
    
    % southeasterly winds, mean: 123 deg from gaussian fit, +-30 deg  
    subplot(223), hold on, box on

    dscatter(bee_dataset.wspd_ms_EWS(ind_SE), plot_var(ind_SE))
    plot_mean_std(bee_dataset.wspd_ms_EWS(ind_SE),plot_var(ind_SE),edges,plot_column)
    
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

    dscatter(bee_dataset.wspd_ms_EWS(ind_rest), plot_var(ind_rest))
    plot_mean_std(bee_dataset.wspd_ms_EWS(ind_rest),plot_var(ind_rest),edges,plot_column)
    
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

if sonde_dT_corr
    
%     figure
%     dscatter(bee_dataset.T_PWS(ind), plot_var(ind))
    
    figure
    txt_pos=0.92;
    edges=-5:2:25;

    % all data
    subplot(221), hold on, box on
    dscatter(bee_dataset.sonde_dT, plot_var)
    
    plot_mean_std(bee_dataset.sonde_dT,plot_var,edges,plot_column)
    
    if plot_column, ylim([0,16]*1e13), else ylim([0.17,1]), end
    xlim([-5.88,24])    
    ylim([0,8]*1e13)
    if plot_column
        ylabel('BrO VCD_{0-4 km} (molec/cm^2)')
    else
         ylabel('VCD_{0-0.6 km} / VCD_{0-4 km}')
    end

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Northerly winds only
    txt_pos=0.92;
    subplot(222), hold on, box on
    ind=(ind_N);
    dscatter(bee_dataset.sonde_dT(ind), plot_var(ind))
    
    plot_mean_std(bee_dataset.sonde_dT(ind),plot_var(ind),edges,plot_column)
    
    if plot_column, ylim([0,16]*1e13), else ylim([0.17,1]), end
    xlim([-5.88,24])    
    ylim([0,8]*1e13)

    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Southeasterly winds only
    subplot(223), hold on, box on
    ind=(ind_SE);
    dscatter(bee_dataset.sonde_dT(ind), plot_var(ind))
    
    plot_mean_std(bee_dataset.sonde_dT(ind),plot_var(ind),edges,plot_column)
    
    if plot_column, ylim([0,16]*1e13), else ylim([0.17,1]), end
    xlim([-5.88,24])    
    ylim([0,8]*1e13)
    xlabel('Sonde dT (°C)')
    if plot_column
        ylabel('BrO VCD_{0-4 km} (molec/cm^2)')
    else
         ylabel('VCD_{0-0.6 km} / VCD_{0-4 km}')
    end

    text(0.95,txt_pos,'Wind: 123° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % all other wind directions
    subplot(224), hold on, box on
    ind=(ind_rest);
    dscatter(bee_dataset.sonde_dT(ind), plot_var(ind))
    
    plot_mean_std(bee_dataset.sonde_dT(ind),plot_var(ind),edges,plot_column)
    
    if plot_column, ylim([0,16]*1e13), else ylim([0.17,1]), end
    xlim([-5.88,24])    
    ylim([0,8]*1e13)
    xlabel('Sonde dT (°C)')

    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    
end

if sonde_dT_o3
    
    plot_var=bee_dataset.o3_surf;
    
    figure
    txt_pos=0.92;
    edges=-5:2:25;

    % all data
    subplot(221), hold on, box on
    ind=~isnan(plot_var);
    dscatter(bee_dataset.sonde_dT(ind), plot_var(ind))
    
    plot_mean_std(bee_dataset.sonde_dT(ind),plot_var(ind),edges,plot_column)
    
    ylim([0,50])
    xlim([-5.88,24])    
    ylabel('surf O_3')

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Northerly winds only
    txt_pos=0.92;
    subplot(222), hold on, box on
    ind=(ind_N & ~isnan(plot_var));
    dscatter(bee_dataset.sonde_dT(ind), plot_var(ind))
    
    plot_mean_std(bee_dataset.sonde_dT(ind),plot_var(ind),edges,plot_column)
    
    ylim([0,50])
    xlim([-5.88,24])    

    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Southeasterly winds only
    subplot(223), hold on, box on
    ind=(ind_SE & ~isnan(plot_var));
    dscatter(bee_dataset.sonde_dT(ind), plot_var(ind))
    
    plot_mean_std(bee_dataset.sonde_dT(ind),plot_var(ind),edges,plot_column)
    
    ylim([0,50])
    xlim([-5.88,24])    
    ylabel('surf O_3')
    xlabel('Sonde dT (°C)')

    text(0.95,txt_pos,'Wind: 123° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % all other wind directions
    subplot(224), hold on, box on
    ind=(ind_rest & ~isnan(plot_var));
    dscatter(bee_dataset.sonde_dT(ind), plot_var(ind))
    
    plot_mean_std(bee_dataset.sonde_dT(ind),plot_var(ind),edges,plot_column)
    
    ylim([0,50])
    xlim([-5.88,24])    
    xlabel('Sonde dT (°C)')

    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    
end

if T_rl_corr
    
%     figure
%     dscatter(bee_dataset.T_PWS(ind), plot_var(ind))
    
    figure
    txt_pos=0.92;

    % all data
    subplot(221), hold on, box on
    dscatter(bee_dataset.T_PWS(ind_t_rl), plot_var(ind_t_rl))
    
    plot_fit_line(bee_dataset.T_PWS(ind_t_rl),plot_var(ind_t_rl),text_size)
    
    ylim([0,8]*1e13)
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
    
    ylim([0,8]*1e13)
    xlim([-45,5])    

    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Southeasterly winds only
    subplot(223), hold on, box on
    ind=(ind_t_rl & ind_SE);
    dscatter(bee_dataset.T_PWS(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.T_PWS(ind),plot_var(ind),text_size)    
    
    ylim([0,8]*1e13)
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
    
    ylim([0,8]*1e13)
    xlim([-45,5])    
    xlabel('PWS T (°C)')

    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    
end


if T_ews_corr
    
%     figure
%     dscatter(bee_dataset.T_PWS(ind), plot_var(ind))
    
    figure
    txt_pos=0.92;

    % all data
    subplot(221), hold on, box on
    dscatter(bee_dataset.T_EWS(ind_t_ews), plot_var(ind_t_ews))
    
    plot_fit_line(bee_dataset.T_EWS(ind_t_ews),plot_var(ind_t_ews),text_size)
    
    ylim([0,8]*1e13)
    xlim([-45,5])    
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Northerly winds only
    txt_pos=0.92;
    subplot(222), hold on, box on
    ind=(ind_t_ews & ind_N);
    dscatter(bee_dataset.T_EWS(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.T_EWS(ind),plot_var(ind),text_size)    
    
    ylim([0,8]*1e13)
    xlim([-45,5])    

    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Southeasterly winds only
    subplot(223), hold on, box on
    ind=(ind_t_ews & ind_SE);
    dscatter(bee_dataset.T_EWS(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.T_EWS(ind),plot_var(ind),text_size)    
    
    ylim([0,8]*1e13)
    xlim([-45,5])    
    xlabel('EWS T (°C)')
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'Wind: 123° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % all other wind directions
    subplot(224), hold on, box on
    ind=(ind_t_ews & ind_rest);
    dscatter(bee_dataset.T_EWS(ind), plot_var(ind))
    
    plot_fit_line(bee_dataset.T_EWS(ind),plot_var(ind),text_size)    
    
    ylim([0,8]*1e13)
    xlim([-45,5])    
    xlabel('EWS T (°C)')

    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    
end

if aod_corr
    
    figure

    subplot(221), hold on
    dscatter(log(bee_dataset.aer_ext),plot_var);
    plot_mean_std(log(bee_dataset.aer_ext),plot_var,1,plot_column)
    ylim([0,8]*1e13)
    xlim([-5,2])
     
    set(gca,'XTick',log([0.01,0.1,1]))
    set(gca,'XMinorTick','on')%log([[0.01:0.01:0.09],[0.1:0.1:0.9],[1:1:2]]))
    tmp=gca;
    tmp.XAxis.MinorTickValues = log([[0.007,0.008,0.009],[0.01:0.01:0.09],[0.1:0.1:0.9],[1:1:7]]);
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    subplot(222), hold on
    dscatter(log(bee_dataset.aer_ext(ind_N)),plot_var(ind_N));
    plot_mean_std(log(bee_dataset.aer_ext(ind_N)),plot_var(ind_N),1,plot_column)
    ylim([0,8]*1e13)
    xlim([-5,2])

    set(gca,'XTick',log([0.01,0.1,1]))
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    subplot(223), hold on
    dscatter(log(bee_dataset.aer_ext(ind_SE)),plot_var(ind_SE));
    plot_mean_std(log(bee_dataset.aer_ext(ind_SE)),plot_var(ind_SE),1,plot_column)
    ylim([0,8]*1e13)
    xlim([-5,2])

    set(gca,'XTick',log([0.01,0.1,1]))
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    subplot(224), hold on
    dscatter(log(bee_dataset.aer_ext(ind_rest)),plot_var(ind_rest));
    plot_mean_std(log(bee_dataset.aer_ext(ind_rest)),plot_var(ind_rest),1,plot_column)
    ylim([0,8]*1e13)
    xlim([-5,2])
   
    set(gca,'XTick',log([0.01,0.1,1]))
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
end

if aod_to_aer
    
    aer_val=bee_dataset.aer_halfmicron;
%     aer_val=bee_dataset.SMPS_100_500;
    figure

    hold on
    ind=~isnan(aer_val);
    dscatter(log(bee_dataset.aer_ext(ind)),aer_val(ind));
    plot_mean_std(log(bee_dataset.aer_ext(ind)),aer_val(ind),1,plot_column)
%     ylim([0,5])
    xlim([-5,2])
     
    set(gca,'XTick',log([0.01,0.1,1]))
    set(gca,'XMinorTick','on')%log([[0.01:0.01:0.09],[0.1:0.1:0.9],[1:1:2]]))
    tmp=gca;
    tmp.XAxis.MinorTickValues = log([[0.007,0.008,0.009],[0.01:0.01:0.09],[0.1:0.1:0.9],[1:1:7]]);
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
    figure

    hold on
    ind=~isnan(aer_val);
    dscatter(aer_val(ind),log(bee_dataset.aer_ext(ind)));
    plot_mean_std(aer_val(ind),log(bee_dataset.aer_ext(ind)),1,plot_column)
%     xlim([0,5])
    ylim([-5,2])
     
    set(gca,'YTick',log([0.01,0.1,1]))
    set(gca,'YMinorTick','on')%log([[0.01:0.01:0.09],[0.1:0.1:0.9],[1:1:2]]))
    tmp=gca;
    tmp.YAxis.MinorTickValues = log([[0.007,0.008,0.009],[0.01:0.01:0.09],[0.1:0.1:0.9],[1:1:7]]);
    set(gca,'YTickLabel',{'10^{-2}','10^{-1}','10^{0}'})
    
end

if aod_wspd
    
    figure
    subplot(221), hold on
    ind=~isnan(bee_dataset.N_SE_rest);
    dscatter(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)))
    plot_mean_std(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)),1)
    xlim([0,15])
    ylim(log([0.014,5]))

    subplot(222), hold on
    ind=bee_dataset.N_SE_rest==1;
    dscatter(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)))
    plot_mean_std(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)),1)
    xlim([0,15])
    ylim(log([0.014,5]))

    subplot(223), hold on
    ind=bee_dataset.N_SE_rest==2;
    dscatter(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)))
    plot_mean_std(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)),1)
    xlim([0,15])
    ylim(log([0.014,5]))

    subplot(224), hold on
    ind=bee_dataset.N_SE_rest==3;
    dscatter(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)))
    plot_mean_std(bee_dataset.wspd_ms(ind),log(bee_dataset.aer_ext(ind)),1)
    xlim([0,15])
    ylim(log([0.014,5]))
    
    

    
end

if aer_wspd
    
    figure
    subplot(221), hold on
    ind=(~isnan(bee_dataset.N_SE_rest) & ~isnan(bee_dataset.aer_halfmicron));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind),1)
    xlim([0,15])
    ylim([0,5])

    subplot(222), hold on
    ind=(bee_dataset.N_SE_rest==1 & ~isnan(bee_dataset.aer_halfmicron));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind),1)
    xlim([0,15])
    ylim([0,5])

    subplot(223), hold on
    ind=(bee_dataset.N_SE_rest==2 & ~isnan(bee_dataset.aer_halfmicron));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind),1)
    xlim([0,15])
    ylim([0,5])

    subplot(224), hold on
    ind=(bee_dataset.N_SE_rest==3 & ~isnan(bee_dataset.aer_halfmicron));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.aer_halfmicron(ind),1)
    xlim([0,15])
    ylim([0,5])
    
end

if smps_wspd
    
    figure
    subplot(221), hold on
    ind=(~isnan(bee_dataset.N_SE_rest) & ~isnan(bee_dataset.SMPS_100_500));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.SMPS_100_500(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.SMPS_100_500(ind),1)
    xlim([0,15])
%     ylim([0,5])

    subplot(222), hold on
    ind=(bee_dataset.N_SE_rest==1 & ~isnan(bee_dataset.SMPS_100_500));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.SMPS_100_500(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.SMPS_100_500(ind),1)
    xlim([0,15])
%     ylim([0,5])

    subplot(223), hold on
    ind=(bee_dataset.N_SE_rest==2 & ~isnan(bee_dataset.SMPS_100_500));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.SMPS_100_500(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.SMPS_100_500(ind),1)
    xlim([0,15])
%     ylim([0,5])

    subplot(224), hold on
    ind=(bee_dataset.N_SE_rest==3 & ~isnan(bee_dataset.SMPS_100_500));
    dscatter(bee_dataset.wspd_ms(ind),bee_dataset.SMPS_100_500(ind))
    plot_mean_std(bee_dataset.wspd_ms(ind),bee_dataset.SMPS_100_500(ind),1)
    xlim([0,15])
%     ylim([0,5])
    
end

if ssa_corr

    figure
    set(gcf, 'Position', [100, 100, 900, 750]);
    fig_ax = tight_subplot(2,2,[0.1,0.07],[0.12,0.07],[0.11,0.05]);
    
    txt_pos=0.92;

    plot_type='hm';
    
    switch plot_type
        case 'sm'
            plot_data=bee_dataset.aer_supermicron;
            xlim_end=1;
            x_label='D_P > 1 \mum (cm^{-3})';
            ind_ssa=(bee_dataset.aer_halfmicron <= 100);
        case 'hm'
            plot_data=bee_dataset.aer_halfmicron;
            xlim_end=6;
            x_label='D_P > 0.5 \mum (cm^{-3})';
            ind_ssa=(bee_dataset.aer_halfmicron <= 100);
    end
    
    % all aer data
    axes(fig_ax(1))
    dscatter(plot_data(ind_ssa), plot_var(ind_ssa)), hold on, box on
    
    plot_fit_line(plot_data(ind_ssa),plot_var(ind_ssa),text_size)
    
    ylim([0,8]*1e13)
    xlim([0,xlim_end])    
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'Wind: \bf{all}', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Northerly winds only
    axes(fig_ax(2))
    ind=(ind_ssa & ind_N);
    dscatter(plot_data(ind), plot_var(ind)), hold on, box on
    
    plot_fit_line(plot_data(ind),plot_var(ind),text_size)    
    
    ylim([0,8]*1e13)
    xlim([0,xlim_end])    

    text(0.95,txt_pos,'Wind: \bf{N}', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % Southeasterly winds only
    axes(fig_ax(3))
    ind=(ind_ssa & ind_SE);
    dscatter(plot_data(ind), plot_var(ind)), hold on, box on
    
    plot_fit_line(plot_data(ind),plot_var(ind),text_size)    
    
    ylim([0,8]*1e13)
    xlim([0,xlim_end])    
    xlabel(x_label)
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    text(0.95,txt_pos,'Wind: \bf{SE}', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    % all other wind directions
    axes(fig_ax(4))
    ind=(ind_ssa & ind_rest);
    dscatter(plot_data(ind), plot_var(ind)), hold on, box on
    
    plot_fit_line(plot_data(ind),plot_var(ind),text_size)    
    
    ylim([0,8]*1e13)
    xlim([0,xlim_end])    
    xlabel(x_label)

    text(0.95,txt_pos,'Wind: \bf{other}', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    
    set(findall(gcf,'-property','FontSize'),'FontSize',17)
    
end

if ssa_corr_oneplot

    figure
    set(gcf, 'Position', [100, 100, 300, 225]);
    fig_ax = tight_subplot(1,1,[0.1,0.07],[0.15,0.09],[0.12,0.07]);
    
    txt_pos=0.92;

    plot_type='hm';
    
    switch plot_type
        case 'sm'
            plot_data=bee_dataset.aer_supermicron;
            xlim_end=1;
            x_label='D_P > 1 \mum (cm^{-3})';
            ind_ssa=(bee_dataset.aer_halfmicron <= 100);
        case 'hm'
            plot_data=bee_dataset.aer_halfmicron;
            xlim_end=6;
            x_label='D_P > 0.5 \mum (cm^{-3})';
            ind_ssa=(bee_dataset.aer_halfmicron <= 100);
    end
    
    % all aer data
    axes(fig_ax(1))
    dscatter(plot_data(ind_ssa), plot_var(ind_ssa)), hold on, box on
    
    plot_fit_line(plot_data(ind_ssa),plot_var(ind_ssa),text_size,'')
    
    ylim([0,12]*1e13)
    xlim([0,xlim_end])    
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel('BrO')
    xlabel('Aer')

%     text(0.95,txt_pos,'Wind: \bf{all}', 'color','k','Units','normalized',...
%         'fontsize',text_size,'HorizontalAlignment','right')
        
    set(findall(gcf,'-property','FontSize'),'FontSize',11)
    
    lines = findobj(gcf,'Type','Line');
    lines(1).LineWidth = 2.0;

end

if smps_corr
    
    figure, hold on
    set(gcf, 'Position', [100, 100, 900, 750]);

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
    
% %     figure, hold on
% %     subplot(221)
% %     dscatter(plot_var(ind_o3),bee_dataset.o3_surf(ind_o3)), hold on
% % %     plot_fit_line(bee_dataset.o3_surf(ind_o3),plot_var(ind_o3),text_size,'right')
% %     xlim([0,7]*1e13)
% %     
% %     subplot(222)
% %     ind=(ind_o3 & ind_N);
% %     dscatter(plot_var(ind),bee_dataset.o3_surf(ind)), hold on
% % %     plot_fit_line(plot_var(ind),bee_dataset.o3_surf(ind),text_size,'right')
% %     xlim([0,7]*1e13)
% % 
% %     subplot(223)
% %     ind=(ind_o3 & ind_SE);
% %     dscatter(plot_var(ind),bee_dataset.o3_surf(ind)), hold on
% % %     plot_fit_line(plot_var(ind),bee_dataset.o3_surf(ind),text_size,'right')
% %     xlim([0,7]*1e13)
% % 
% %     subplot(224)
% %     ind=(ind_o3 & ind_rest);
% %     dscatter(plot_var(ind),bee_dataset.o3_surf(ind)), hold on
% % %     plot_fit_line(bee_dataset.o3_surf(ind),plot_var(ind),text_size,'right')
% %     xlim([0,7]*1e13)
    
    
    
    % box plot
    figure
    nbars=3;
    a=0.2;
    
    data=bee_dataset.bro_col;
    
    o3_var=NaN(size(bee_dataset.o3_surf));
    o3_var(bee_dataset.o3_surf<10)=1;
    o3_var(bee_dataset.o3_surf>=10 & bee_dataset.o3_surf<15)=2;
    o3_var(bee_dataset.o3_surf>=15 & bee_dataset.o3_surf<25)=3;
    o3_var(bee_dataset.o3_surf>=25)=4;
    
    box_group=o3_var*10+bee_dataset.N_SE_rest;

    % get position of each bar
    box_pos=sort([1:a:1+a*(nbars-1),2:a:2+a*(nbars-1),3:a:3+a*(nbars-1),4:a:4+a*(nbars-1)]);

    % set tick position as the mean of each group
    tick_pos=[];
    for i=1:4
        tick_pos=[tick_pos,mean(box_pos((nbars)*(i-1)+1:(nbars)*i))];
    end

    c_list=winter(3);
    % do box plot
    boxplot(data,box_group,'positions',box_pos,'colors',c_list,...
            'jitter',0.5,'symbol','.')
    hold on, box on
        
    % plot mean as well
    group_ind=unique(box_group(~isnan(box_group)));
    c_ind=repmat(1:nbars,1,length(group_ind)/nbars);
    for i=1:length(group_ind)
        plot(box_pos(i), nanmean(data(box_group==group_ind(i))),'x','color',c_list(c_ind(i),:))
    end
        
    % some formatting
    set(findobj(gca,'type','line'),'linew',1.5)
    set(findall(gca,'tag','Outliers'),'MarkerSize',4);
    
    set(gca,'xtick',tick_pos)
    set(gca,'xticklabel',{'<10','10-15','15-25','>25'})
    
    xlim([box_pos(1)-0.15,box_pos(end)+0.15])
    
    set(gca, 'YGrid', 'on')

    % legend for one figure only -- position set later
    ll=legend(flipud(findall(gca,'Tag','Box')),...
           {'N','SE','other'},...
           'Orientation','horizontal','location','northeast');
    
    
       
       
% %     % bar plot of BrO in ozone bins (ignores NaN)
% %     bro_limit=mean(bee_dataset.bro_col);
% % %     bro_sort=sort(bee_dataset.bro_col);
% % %     bro_limit=bro_sort(floor(length(bee_dataset.bro_col)*0.75)+1);
% %     
% %     o3_0_10=(bee_dataset.o3_surf<10 & bee_dataset.bro_col>=bro_limit);
% %     o3_10_20=(bee_dataset.o3_surf>=10 & bee_dataset.o3_surf<20 & ...
% %               bee_dataset.bro_col>=bro_limit);
% %     o3_20_30=(bee_dataset.o3_surf>=20 & bee_dataset.o3_surf<30 & ...
% %               bee_dataset.bro_col>=bro_limit);
% %     o3_30_up=(bee_dataset.o3_surf>=30 & bee_dataset.bro_col>=bro_limit);
% %     
% %     by_wind=[ [sum(o3_0_10 & ind_N),...
% %                 sum(o3_0_10 & ind_SE),...
% %                 sum(o3_0_10 & ind_rest) ];...
% %               [sum(o3_10_20 & ind_N),...
% %                 sum(o3_10_20 & ind_SE),...
% %                 sum(o3_10_20 & ind_rest) ];...
% %               [sum(o3_20_30 & ind_N),...
% %                 sum(o3_20_30 & ind_SE),...
% %                 sum(o3_20_30 & ind_rest) ];...
% %               [sum(o3_30_up & ind_N),...
% %                 sum(o3_30_up & ind_SE),...
% %                 sum(o3_30_up & ind_rest) ] ];
% %             
% %     figure
% % %     bar(by_wind./sum(sum(by_wind)),'stacked');
% %     bar(by_wind./sum(sum(by_wind)));
% %     legend('354° \pm 30°','123° \pm 30°','Other','location','north')
% %     xlabel('Surface O_3 (ppbv)')
% %     ylabel('p for BrO VCD_{0-4 km} > mean')
% %     xlim([0.5,4.5])
% %     
% %     c = {'0-10','10-20','20-30','>30'};
% %     set(gca,'xticklabel',c)
    
    
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

if SI_corr
    
    edges=[0:0.5:11]*1e13;
    
    txt_pos=0.92;
         
    plot_x=bee_dataset.FYSI_3day;
    
    si_x_lim=log([50,8.1e5]);
    
    % all winds
    figure
    subplot(221), hold on, box on
    
    ind_ok=~isnan(plot_x);
    dscatter(plot_x(ind_ok), plot_var(ind_ok))
%     plot_mean_std(plot_x(ind_ok),plot_var(ind_ok),edges)

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    ylim([0,12]*1e13)
    xlim(si_x_lim)    
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')
    
    % northerly winds, mean: 354 deg from gaussian fit, +-30 deg
    subplot(222), hold on, box on

    ind=(ind_ok & ind_N);
    dscatter(plot_x(ind), plot_var(ind))
%     plot_mean_std(plot_x(ind),plot_var(ind),edges)
%     plot_fit_line(plot_x(ind),plot_var(ind),12)
    
    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    ylim([0,12]*1e13)
    xlim(si_x_lim)    
    
    % southeasterly winds, mean: 123 deg from gaussian fit, +-30 deg  
    subplot(223), hold on, box on

    ind=(ind_ok & ind_SE);
    dscatter(plot_x(ind), plot_var(ind))
%     plot_mean_std(plot_x(ind),plot_var(ind),edges)
    
    text(0.95,txt_pos,'Wind: 123° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    ylim([0,12]*1e13)
    xlim(si_x_lim)    
    xlabel('FYSI contact (s)')
    ylabel('BrO VCD_{0-4 km} (molec/cm^2)')

    % everything else
    subplot(224), hold on, box on

    ind=(ind_ok & ind_rest);
    dscatter(plot_x(ind), plot_var(ind))
%     plot_mean_std(plot_x(ind),plot_var(ind),edges)
    
    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    ylim([0,12]*1e13)
    xlim(si_x_lim)    
    xlabel('FYSI contact (s)')
    
    
end

if SI_corr_flip
    
    edges=[0:0.5:11]*1e13;
    
    txt_pos=0.92;
         
    plot_x=plot_var;
    plot_var=bee_dataset.FYSI_3day;
    
    % all winds
    figure
    subplot(221), hold on, box on

    ind_ok=~isnan(plot_var);
    dscatter(plot_x(ind_ok), plot_var(ind_ok))
    plot_mean_std(plot_x(ind_ok),plot_var(ind_ok),edges)

    text(0.95,txt_pos,'Wind: 0-360°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    ylim([0,11]*1e13)
    xlim([0,12]*1e13)   
    ylabel('FYSI contact (s)')
    
    % northerly winds, mean: 354 deg from gaussian fit, +-30 deg
    subplot(222), hold on, box on

    ind=(ind_ok & ind_N);
    dscatter(plot_x(ind), plot_var(ind))
    plot_mean_std(plot_x(ind),plot_var(ind),edges)
    
    text(0.95,txt_pos,'Wind: 354° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    ylim([0,11]*1e13)
    xlim([0,12]*1e13)   
    
    % southeasterly winds, mean: 123 deg from gaussian fit, +-30 deg  
    subplot(223), hold on, box on

    ind=(ind_ok & ind_SE);
    dscatter(plot_x(ind), plot_var(ind))
    plot_mean_std(plot_x(ind),plot_var(ind),edges)
    
    text(0.95,txt_pos,'Wind: 123° \pm 30°', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    ylim([0,11]*1e13)
    xlim([0,12]*1e13)    
    xlabel('BrO VCD_{0-4 km} (molec/cm^2)')
    ylabel('FYSI contact (s)')

    % everything else
    subplot(224), hold on, box on

    ind=(ind_ok & ind_rest);
    dscatter(plot_x(ind), plot_var(ind))
    plot_mean_std(plot_x(ind),plot_var(ind),edges)
    
    text(0.95,txt_pos,'Wind: other', 'color','k','Units','normalized',...
        'fontsize',text_size,'HorizontalAlignment','right')
    ylim([0,11]*1e13)
    xlim([0,12]*1e13)    
    xlabel('BrO VCD_{0-4 km} (molec/cm^2)')
    
    
end

if SI_to_SI
    
    ind_ok=(~isnan(bee_dataset.FYSI_3day) & ~isnan(bee_dataset.MYSI_3day));

    figure
    loglog(bee_dataset.FYSI_3day(ind_ok)+bee_dataset.MYSI_3day(ind_ok),...
         bee_dataset.FYSI_3day(ind_ok),'b.')
    xlabel('Total ice contact (s)')
    ylabel('FYI contact (s)')
    
    figure

    subplot(221), box on
    ind=ind_ok;
    loglog(bee_dataset.FYSI_3day(ind),bee_dataset.MYSI_3day(ind),'b.')
    ylabel('MYSI contact (s)')

    subplot(222), box on
    ind=(ind_ok & ind_N);
    loglog(bee_dataset.FYSI_3day(ind),bee_dataset.MYSI_3day(ind),'b.')

    subplot(223), box on
    ind=(ind_ok & ind_SE);
    loglog(bee_dataset.FYSI_3day(ind),bee_dataset.MYSI_3day(ind),'b.')
    xlabel('FYSI contact (s)')
    ylabel('MYSI contact (s)')

    subplot(224), box on
    ind=(ind_ok & ind_rest);
    loglog(bee_dataset.FYSI_3day(ind),bee_dataset.MYSI_3day(ind),'b.')
    xlabel('FYSI contact (s)')
    
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

if windrose
   
    % PWS wind data
    load('/home/kristof/work/weather_stations/ridge_lab/PWS_all.mat');
    data((month(data.DateTime)>5 | month(data.DateTime)<3),:)=[];

    WindRose(data.WindDir,data.WindSpd,'AngleNorth',0,'AngleEast',90,...
             'FreqLabelAngle',45,'vWinds',[0 3 6 9 12 15]);
    
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

function plot_mean_std(xx,yy,edges,plot_column)

    if nargin==3, plot_column=1; end
    
    tmp_mean=[];
    tmp_std=[];
    
    if length(edges)==1
        % use deciles
        edges=prctile(xx,0:10:100);
    end

    for i=1:length(edges)-1
        tmp_mean=[tmp_mean,nanmean(yy(xx>=edges(i) & xx < edges(i+1)))];
        tmp_std=[tmp_std,nanstd(yy(xx>=edges(i) & xx < edges(i+1)))];
    end    
    
    if ~plot_column, plot([edges(1),edges(end)],[0.83,0.83],'b-','linewidth',1), end

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
        
    else
        % no text
    end
    
end





