% Principal component analysis for BEE data

%% setup

weather_only=0; % use weather data only (adds OPC and SMPS if set to false)
weather_comp=0; % weather only, but remove 2015 to compare to aer data

plot_pca=0; % all components wth var>1, on separate plots by wind speed
plot_pca_poster=1; % first 2 PC for N and SE winds on one plot

% load dataset
load('/home/kristof/work/BEEs/BEE_dataset_all.mat')

% variables included in PCA
if weather_only
    vars={'BrO','Ext','V','T','dT_{200}','dT_{600}','PBLH','P','\DeltaP'};
else
    vars={'BrO','Ext','V','T','dT_{200}','dT_{600}','PBLH','P','\DeltaP','Aer','FYSI','MYSI'};
%     vars={'BrO','Ext','V','T','dT_{200}','dT_{600}','P','\DeltaP','Aer','SMPS'};
end

% wind direction key
wdir_key={'all','\bf{N}','\bf{SE}','other'};
wdir_key2={'Wind: all','Wind: \bf{N}','Wind: \bf{SE}','Wind: other'};

%% loopover each wind direction

if plot_pca_poster
    figure(99)
    set(gcf, 'Position', [100, 100, 1100, 310]);
    fig_ax = tight_subplot(2,2,[0.15,0.04],[0.11,0.1],[0.05,0.01]);
end

for wdir=0:3
    
    % select data by wind direction
    if wdir==0
        ind=true(size(bee_dataset.bro_col));
    else
        ind=bee_dataset.N_SE_rest==wdir;
    end
    
    
    %% create data arrays for pca
        
    % weather conditions only
    data_in_weather=[bee_dataset.bro_col(ind),...
                     bee_dataset.aer_ext(ind),...
                     bee_dataset.wspd_ms(ind),...
                     bee_dataset.T_PWS(ind),...
                     bee_dataset.sonde_dT(ind),...
                     bee_dataset.sonde_T_200(ind)-bee_dataset.sonde_T_0(ind),...
                     bee_dataset.bl_height(ind),...
                     bee_dataset.P_PWS(ind),...
                     bee_dataset.P_PWS_tend(ind),...
                     ];

    % add aerosol data
    data_in_aer=[bee_dataset.bro_col(ind),...
                 bee_dataset.aer_ext(ind),...
                 bee_dataset.wspd_ms(ind),...
                 bee_dataset.T_PWS(ind),...
                 bee_dataset.sonde_dT(ind),...
                 bee_dataset.sonde_T_200(ind)-bee_dataset.sonde_T_0(ind),...
                 bee_dataset.bl_height(ind),...
                 bee_dataset.P_PWS(ind),...
                 bee_dataset.P_PWS_tend(ind),...
                 bee_dataset.aer_halfmicron(ind),...
                 bee_dataset.FYSI_3day(ind),...
                 bee_dataset.MYSI_3day(ind),...
                 ];
%     data_in_aer=[bee_dataset.bro_col(ind),...
%                  bee_dataset.aer_ext(ind),...
%                  bee_dataset.wspd_ms(ind),...
%                  bee_dataset.T_PWS(ind),...
%                  bee_dataset.sonde_dT(ind),...
%                  bee_dataset.sonde_T_200(ind)-bee_dataset.sonde_T_0(ind),...
%                  bee_dataset.P_PWS(ind),...
%                  bee_dataset.P_PWS_tend(ind),...
%                  bee_dataset.aer_halfmicron(ind),...
%                  bee_dataset.SMPS_100_500(ind),...
%                  ];

    data_times=bee_dataset.times(ind);
             
    if weather_only
        
        % weather data only
        data_in=data_in_weather;
        
        if weather_comp
            % weather data only, but remove 2015 to compare to aerosol results
            data_in=data_in_aer;
            data_in(data_times.Year==2015,:)=[];
            data_in(:,end-1:end)=[];
        end
        
    else
        
        % weather data plus OPC and SMPS
        data_in=data_in_aer;
        
    end
    
    %% do PCA
    
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
    [coeffs,~,pc_var,~,tot_var] = pca(data_in,'Centered',false);

    % number of components with variance>x
    n=sum(pc_var>1);

    %% plot PCA results
    if plot_pca
    
        figure, hold on

        for i=1:n

            subplot(n,1,i)

            % loading greater than cutoff, for coloring bars
            ind=abs(coeffs(:,i))>0.25;
            
            % flip so BrO is consistently positive (easier to interpret
            % results -- sign is arbitrary anyway)
            if coeffs(1,i)<0, coeffs(:,i)=-coeffs(:,i); end

            if sum(ind)==0
                % no loadings greater than cutoff
                bar(coeffs(:,i));
            else

                % plot loadings above and below the cutoff separately
                x=1:length(vars);

                y1=coeffs(:,i);
                y1(ind)=0; % < cutoff
                y2=coeffs(:,i);
                y2(~ind)=0; % > cutoff

                b1=bar(x,y1); hold on
                b2=bar(x,y2);

                b1.FaceAlpha = 0.5;


            end

            ylim([-0.8,0.8])
            xlim([0.5,length(vars)+0.5])

            ylabel(['PC' num2str(i)])

            % add percent of variance explained to plot
            text(0.01,0.8,['Var: ' num2str(round(tot_var(i),1)) '%'],...
                'color','k','Units','normalized','HorizontalAlignment','Left')

            % use labels for last subplot only
            if i==n
                set(gca,'XTickLabel',vars);
            else
                set(gca,'XTickLabel',[]);

                if i==1, title(wdir_key2{wdir+1},'fontweight','normal'); end
            end

        end
    end
    
    if plot_pca_poster
        
        if wdir==0 || wdir==3
            continue
        end
            
        figure(99)
        
        for i=1:2

            plotind=2*i+wdir-2;
            axes(fig_ax(plotind))

            % loading greater than cutoff, for coloring bars
            ind=abs(coeffs(:,i))>0.25;
            
            % flip so BrO is consistently positive (easier to interpret
            % results -- sign is arbitrary anyway)
            if coeffs(1,i)<0, coeffs(:,i)=-coeffs(:,i); end

            % plot loadings above and below the cutoff separately
            x=1:length(vars);

            y1=coeffs(:,i);
            y1(ind)=0; % < cutoff
            y2=coeffs(:,i);
            y2(~ind)=0; % > cutoff

            b1=bar(x,y1); hold on
            b2=bar(x,y2);

            b1.FaceAlpha = 0.5;

            ylim([-0.6,0.6])
            xlim([0.5,length(vars)+0.5])

            if any(plotind==[1,3]), ylabel(['PC' num2str(i)]); end

%             % add percent of variance explained to plot
%             text(0.01,0.8,['Var: ' num2str(round(tot_var(i),1)) '%'],...
%                 'color','k','Units','normalized','HorizontalAlignment','Left')

            % use labels for last subplot only
            if i==2
                set(gca,'XTickLabel',vars,'FontSize',10);
            else
%                 set(gca,'XTickLabel',[],'FontSize',11);
                set(gca,'XTickLabel',vars,'FontSize',10);
                title(wdir_key2{wdir+1},'fontweight','normal','FontSize',17)
            end

            
        end

    end
    
end
