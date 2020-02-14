% plot yearly results form MAX-DOAS profile retrievals

plot_monthly=1; % 1 for monthly box plots, 0 for yearly box plots

plot_bro=0;
plot_aod=0;
plot_dT=0;
plot_T_10=0;
plot_ws_dT=0;
plot_BL=0;
plot_ssa=0;
plot_fysi=0;
plot_o3=0;

% not monthly/yearly box plots:
plot_dT_prob=1;

if any([plot_dT_prob]==1), plot_monthly=3; end

%% load/filter BEE dataset
load('/home/kristof/work/BEEs/BEE_dataset_all.mat')

% exclude 2015 (intermittent measurements), and April 2017 (4 days only)
filt_ind=find(bee_dataset.times.Year==2015 | ...
              (bee_dataset.times.Year==2017 & bee_dataset.times.Month==4)...
              );

bee_dataset(filt_ind,:)=[];


%% setup
c={'r','g','b','m','c'};
m_names={'March','April','May'};

plot_x=bee_dataset.times;

if plot_bro
%     tmp=(bee_dataset.wspd_ms_EWS<2 & bee_dataset.o3_surf<15 & bee_dataset.sonde_dT<7);
%     plot_y=bee_dataset.bro_col(tmp);
%     plot_x=plot_x(tmp);
    
    plot_y=bee_dataset.bro_col;
    y_label='BrO VCD (molec/cm^2)';
    yscale='linear';
    ymin=-5e12;
end

if plot_aod
    plot_y=bee_dataset.aer_ext;
    y_label='AOD';
    yscale='log';
    ymin=0.01;
end

if plot_dT
    t_sonde=[];
    for yr=unique(bee_dataset.times.Year)'
        t_sonde=[t_sonde; get_sonde_PT(yr)];
    end
    
    t_sonde(t_sonde.date.Month<3 | t_sonde.date.Month>5,:)=[];
    
    plot_y=t_sonde.dT;
    plot_x=t_sonde.date;

    y_label='Radiosonde T_{600m} - T_{10m} (^{\circ}C)';
    yscale='linear';
    ymin=-6;
end

if plot_ws_dT
    
    % PWS data
    load('/home/kristof/work/weather_stations/ridge_lab/PWS_all.mat');
    data((month(data.DateTime)>5 | month(data.DateTime)<3),:)=[];
    data(year(data.DateTime)==2015,:)=[];
    pws_data=data;

    % EWS data
    load('/home/kristof/work/weather_stations/Eureka/EWS_PTU_and_weather_complete.mat')
    data(year(data.DateTime)==2015,:)=[];
    ews_data=data;    
    
    tmp=find_coincident_mean(ews_data.DateTime,pws_data.DateTime,pws_data.TempC,30);
    
    plot_y=tmp-ews_data.TempC;
    plot_x=ews_data.DateTime;

    y_label='PWS T - EWS T (^{\circ}C)';
    yscale='linear';
    ymin=-8;
    
end

if plot_T_10
    
    load('/home/kristof/work/weather_stations/Eureka/EWS_PTU_and_weather_complete.mat')
    data(data.Year<2016,:)=[];
    
    plot_y=data.TempC;
    plot_x=data.DateTime;

    y_label='T_{10m} (^{\circ}C)';
    yscale='linear';
    ymin=-50;
end

if plot_BL
    BL_out = get_BL_height( unique(bee_dataset.times.Year)' );
    BL_out(BL_out.DateTime.Month<3 | BL_out.DateTime.Month>5,:)=[];
    
    plot_y=BL_out.BL_height_m;
    plot_x=BL_out.DateTime;

    y_label='Radiosonde T_{600m} - T_{10m} (^{\circ}C)';
    yscale='linear';
    ymin=-6;
end

if plot_ssa
    plot_y=bee_dataset.aer_halfmicron;
    y_label='D_p > 0.5 \mum';
    yscale='linear';
    ymin=0;
end

if plot_fysi
    plot_y=bee_dataset.FYSI_3day;
    y_label='3 day FYSI contact (s m^2)';
    yscale='linear';
    ymin=0;
end

if plot_o3
    plot_y=bee_dataset.o3_surf;
    y_label='Surface ozone (ppbv)';
    yscale='linear';
    ymin=0;
end

%% plot monthly stats
if plot_monthly==1
    
    figure
    set(gcf, 'Position', [100, 100, 1000, 500]);

    for mm=3:5
        subplot(1,3,mm-2)
        ind=find(plot_x.Month==mm);

        yrs=plot_x.Year(ind);

        boxplot(plot_y(ind),yrs,'positions',unique(yrs)'-2015,'jitter',0.5), hold on
        set(findall(gca,'tag','Outliers'),'MarkerSize',4);

        if mm==3, ylabel(y_label); end

        for yr=unique(yrs)'
            ind2=yrs==yr;
            plot(yr-2015,nanmean(plot_y(ind(ind2))),'kx','linewidth',1.2)
        end

        title(m_names(mm-2))
        ylim([ymin,max(plot_y)*1.03])
        set(gca, 'YScale', yscale)
        set(gca, 'YGrid', 'on')

    end

    % figure, hold on
    % for yr=2015:2019
    %     
    %     ind=find(plot_x.Year==yr);
    %     plot(plot_x(ind),plot_y(ind),'o','color',c{yr-2014})
    %     
    % end
    
elseif plot_monthly==0
   
    figure
    set(gcf, 'Position', [100, 100, 1200, 650]);
    fig_ax = tight_subplot(4,1,[0.07,0.04],[0.09,0.08],[0.09,0.04]);
    
    for yr=2016:2019
    
        axes(fig_ax(yr-2015))
        
        ind=(plot_x.Year==yr);
        plot(plot_x(ind),plot_y(ind),'rx','markerfacecolor','r','markersize',6)
        
        grid on
        ylim([ymin,max(plot_y)*1.03])
        xlim([datenum(yr,3,1),datenum(yr,6,2)])
        
        set(gca, 'XTick', datenum([[yr,3,1];[yr,3,15];[yr,4,1];[yr,4,15];...
                                   [yr,5,1];[yr,5,15];[yr,6,1]]))
        
        if yr~=2019, set(gca,'xticklabel',[]), end
        set(gca, 'YScale', yscale)
        
        if yr==2017
            ylb=ylabel(y_label); 
            ylb.Position(1)=ylb.Position(1)-3;
            ylb.Position(2)=ylb.Position(2)-max(plot_y);
        end
        
        text(0.9,0.8,num2str(yr), 'color','k','Units','normalized',...
        'fontsize',12,'fontweight','bold')
        
        set(findall(gcf,'-property','FontSize'),'FontSize',17)

    end
    
else
    
    if plot_dT_prob
        % get change of probability of inversion strength <2C when
        % comparing winds above and below 8 m/s
        
        dT_stats=table;
        
        ind_pws1=bee_dataset.wspd_ms<8;
        ind_pws2=bee_dataset.wspd_ms>=8;
        ind_ews1=bee_dataset.wspd_ms_EWS<8;
        ind_ews2=bee_dataset.wspd_ms_EWS>=8;

        bins_edges=-4:6:20;
        
        
        for i=0:3
            
            if i==0
                tmp=true(size(bee_dataset.N_SE_rest));
            else
                tmp=bee_dataset.N_SE_rest==i;
            end
            
            pws1=histcounts(bee_dataset.sonde_dT(ind_pws1 & tmp),bins_edges,'normalization','probability');
            pws2=histcounts(bee_dataset.sonde_dT(ind_pws2 & tmp),bins_edges,'normalization','probability');

            ews1=histcounts(bee_dataset.sonde_dT(ind_ews1 & tmp),bins_edges,'normalization','probability');
            ews2=histcounts(bee_dataset.sonde_dT(ind_ews2 & tmp),bins_edges,'normalization','probability');
            
            pws_up=((pws2(1)/pws1(1))-1)*100;
            ews_up=((ews2(1)/ews1(1))-1)*100;
            
            disp(['PWS and EWS p increase (%), i=' num2str(i)])
            disp([pws_up,ews_up])
            
        end
    end
    
end





