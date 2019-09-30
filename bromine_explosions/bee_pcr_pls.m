% Principal component analysis for BEE data

%% setup

weather_only=0; % use weather data only (adds OPC and SMPS id set to false)
weather_comp=0; % weather only, but remove 2015 to compare to aer data

% load dataset
load('/home/kristof/work/BEEs/BEE_dataset_all.mat')

% variables included in PCA
if weather_only
    vars={'BrO','EXT','V','T','dT_{200}','dT_{600}','P','\DeltaP'};
else
    vars={'BrO','EXT','V','T','dT_{200}','dT_{600}','P','\DeltaP','OPC'};
%     vars={'BrO','EXT','V','T','dT_{200}','dT_{600}','P','\DeltaP','OPC','SMPS'};
end

% wind direction key
wdir_key={'all','N','SE','other'};

%% loopover each wind direction

for wdir=0:3
    
    % select data by wind direction
    if wdir==0
        ind=true(size(bee_dataset.bro_col));
    else
        ind=bee_dataset.N_SE_rest==wdir;
    end
    
    
    %% create data arrays
    
    % weather conditions only
    data_in_weather=[bee_dataset.bro_col(ind),...
                     bee_dataset.aer_ext(ind),...
                     bee_dataset.wspd_ms(ind),...
                     bee_dataset.T_PWS(ind),...
                     bee_dataset.sonde_dT(ind),...
                     bee_dataset.sonde_T_200(ind)-bee_dataset.sonde_T_0(ind),...
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
                 bee_dataset.P_PWS(ind),...
                 bee_dataset.P_PWS_tend(ind),...
                 bee_dataset.OPC_supermicron(ind),...
                 ];
%     data_in_aer=[bee_dataset.bro_col(ind),...
%                  bee_dataset.aer_ext(ind),...
%                  bee_dataset.wspd_ms(ind),...
%                  bee_dataset.T_PWS(ind),...
%                  bee_dataset.sonde_dT(ind),...
%                  bee_dataset.sonde_T_200(ind)-bee_dataset.sonde_T_0(ind),...
%                  bee_dataset.P_PWS(ind),...
%                  bee_dataset.P_PWS_tend(ind),...
%                  bee_dataset.OPC_supermicron(ind),...
%                  bee_dataset.SMPS_100_500(ind),...
%                  ];
                 
    data_times=bee_dataset.times(ind);
             
    if weather_only
        
        % weather data only
        data_in=data_in_weather;
        
        % ptom data for prediction
        data_ptom=[bee_dataset.wspd_ms(ind),...
                   bee_dataset.ptom_T_600(ind),...
                   bee_dataset.ptom_T_600(ind)-bee_dataset.ptom_T_0(ind),...
                   bee_dataset.ptom_T_200(ind)-bee_dataset.ptom_T_0(ind),...
                   bee_dataset.ptom_P_600(ind),...
                   bee_dataset.ptom_P_600_tend(ind),...
                   ];
        
        if weather_comp
            % weather data only, but remove 2015 to compare to aerosol results
            data_in(data_times.Year==2015,:)=[];
            data_ptom(data_times.Year==2015,:)=[];
            data_times(data_times.Year==2015,:)=[];
        end
        
    else
        
        % weather data plus OPC and SMPS
        data_in=data_in_aer;
        
        % ptom data for prediction
        data_ptom=[bee_dataset.wspd_ms(ind),...
                   bee_dataset.ptom_T_600(ind),...
                   bee_dataset.ptom_T_600(ind)-bee_dataset.ptom_T_0(ind),...
                   bee_dataset.ptom_T_200(ind)-bee_dataset.ptom_T_0(ind),...
                   bee_dataset.ptom_P_600(ind),...
                   bee_dataset.ptom_P_600_tend(ind),...
                   bee_dataset.ptom_supermicron(ind),...
                   ];
%         data_ptom=[bee_dataset.wspd_ms(ind),...
%                    bee_dataset.ptom_T_600(ind),...
%                    bee_dataset.ptom_T_600(ind)-bee_dataset.ptom_T_0(ind),...
%                    bee_dataset.ptom_T_200(ind)-bee_dataset.ptom_T_0(ind),...
%                    bee_dataset.ptom_P_600(ind),...
%                    bee_dataset.ptom_P_600_tend(ind),...
%                    bee_dataset.ptom_supermicron(ind),...
%                    bee_dataset.ptom_submicron(ind),...
%                    ];
        
    end
    
    %% predict BrO using partial least squares regression
    
    % remove NaN's for pls regression
    remove_ind=sum(isnan(data_in),2)>0;
    data_tmp=data_in(~remove_ind,:);
    
    xx=data_tmp(:,3:end); % use everything except BrO and extinction to predict
    yy=data_tmp(:,1); % BrO values
    
    % get PLs fit
    [Xloadings,Yloadings,Xscores,Yscores,betaPLS] = plsregress(xx,yy,5);
    
    % get predicted BrO
    yfitPLS = [ones(size(xx,1),1) xx]*betaPLS;
    
    % plot and calculate R^2
%     figure
%     plot([-1,12]*1e13,[-1,12]*1e13,'k--'), hold on
%     dscatter(yy,yfitPLS)
%     title(['PLS; wind direction: ' wdir_key{wdir+1}])
% 
%     xlim([-1,12]*1e13)
%     ylim([-1,12]*1e13)
    
    TSS = sum((yy-mean(yy)).^2);
    RSS_PLS = sum((yy-yfitPLS).^2);
    rsquaredPLS = 1 - RSS_PLS/TSS;
    
    %% predict BrO from pTOMCAT met and aer fields using PLS results
    
    yfitPLS = [ones(size(data_ptom,1),1) data_ptom]*betaPLS;
    
%     figure
%     plot([-1,12]*1e13,[-1,12]*1e13,'k--'), hold on
%     dscatter(yy,yfitPLS(~remove_ind))
%     title(['PLS; wind direction: ' wdir_key{wdir+1}])
% 
%     xlim([-1,12]*1e13)
%     ylim([-1,12]*1e13)
    
    RSS_PLS = sum((yy-yfitPLS(~remove_ind)).^2);
    rsquaredPLS_new = 1 - RSS_PLS/TSS;
    
    %% predict BrO using principal component regression

    % center (mean=0) and standardize (std=1) (setting 'VariableWeights','variance'
    % in pca doesn't actually standardize the data for some reason)
    data_mean=nanmean(data_in);
    data_std=nanstd(data_in);
    
    for i=1:size(data_in,2)
        data_in(:,i)=data_in(:,i)-data_mean(i);
        data_in(:,i)=data_in(:,i)/data_std(i);
    end
    
    xx=data_in(:,3:end); % use everything except BrO and extinction to predict
    yy=data_in(:,1); % BrO values

    [PCALoadings,PCAScores,PCAVar] = pca(xx,'centered',false);
    betaPCR = regress(yy, PCAScores(:,1:5));  
    
    betaPCR = PCALoadings(:,1:5)*betaPCR;
    betaPCR = [nanmean(xx)*betaPCR; betaPCR];
    yfitPCR = [ones(size(xx,1),1) xx]*betaPCR;
    
%     figure
%     plot([-1,12]*1e13,[-1,12]*1e13,'k--'), hold on
%     dscatter((yy(~isnan(yfitPCR))*data_std(1))+data_mean(1),...
%              (yfitPCR(~isnan(yfitPCR))*data_std(1))+data_mean(1))
%     title(['PCR; wind direction: ' wdir_key{wdir+1}])
%     
%     xlim([-1,12]*1e13)
%     ylim([-1,12]*1e13)

    TSS = sum((yy).^2);
    RSS_PCR = nansum((yy-yfitPCR).^2);
    rsquaredPCR = 1 - RSS_PCR/TSS;
    
%     disp(['R^2: PLS: ' num2str(rsquaredPLS) ' PCR: ' num2str(rsquaredPCR)])

    %% predict BrO from pTOMCAT met and aer fields using PCR results

    for i=1:size(data_ptom,2)
        data_ptom(:,i)=data_ptom(:,i)-nanmean(data_ptom(:,i));
        data_ptom(:,i)=data_ptom(:,i)/nanstd(data_ptom(:,i));
    end
    
    yfitPCR = [ones(size(data_ptom,1),1) data_ptom]*betaPCR;
    
    figure
    plot([-1,12]*1e13,[-1,12]*1e13,'k--'), hold on
    dscatter((yy(~isnan(yfitPCR))*data_std(1))+data_mean(1),...
             (yfitPCR(~isnan(yfitPCR))*data_std(1))+data_mean(1))
    title(['PCR; wind direction: ' wdir_key{wdir+1}])
    
    xlim([-1,12]*1e13)
    ylim([-1,12]*1e13)

    TSS = sum((yy).^2);
    RSS_PCR = nansum((yy-yfitPCR).^2);
    rsquaredPCR_new = 1 - RSS_PCR/TSS;
    
    disp(['PLS: ' num2str(rsquaredPLS) 'PLS ptom: ' num2str(rsquaredPLS_new) ...
          ' PCR: ' num2str(rsquaredPCR) ' PCR ptom: ' num2str(rsquaredPCR_new)])
      
    if wdir==0
        bro_meas=(yy*data_std(1))+data_mean(1);
        pca_bro=(yfitPCR*data_std(1))+data_mean(1);
        pca_times=data_times;
        
        figure
        plot(bee_dataset.times,bee_dataset.bro_col_ptom,'ro'), hold on
        plot(data_times,pca_bro,'bo')

        plot(bee_dataset.times,bee_dataset.bro_col,'ks')
        
        save('/home/kristof/work/BEEs/ptom_pcr.mat','pca_times','pca_bro');
        
    end
    
end


