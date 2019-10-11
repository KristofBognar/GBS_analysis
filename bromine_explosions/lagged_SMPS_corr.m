% function [ output_args ] = lagged_SMPS_corr()
%LAGGED_SMPS_CORR Summary of this function goes here
%   Detailed explanation goes here

%% setup

redo_lagged=0;

% load BrO dataset
load('/home/kristof/work/BEEs/BEE_dataset_all.mat')

bee_dataset(year(bee_dataset.times)==2015,:)=[];

% get duration of each profile
prof_len=match_prof_length(bee_dataset.times);
prof_len=prof_len/2; % +- minutes around mean time

% load SMPS data
load('/home/kristof/work/SMPS/smps+opc+aps/smps_size_dist_all.mat');

% integrate Dp data
% 50-500 nm: 23-end
% 100-500 nm: 33-end
smps_sum=sum(smps_data(:,33:end),2);

ind=find(smps_sum==0 | isnan(smps_sum));
smps_data(ind,:)=[];
smps_time(ind)=[];
smps_tot_data(ind)=[];
smps_sum(ind)=[];


%% get lagged aerosol concentrations

if redo_lagged
    
    lags=duration([[1:72]' zeros(72,2)]);

    smps_mean=NaN(length(bee_dataset.times),length(lags));
    R2_array=NaN(length(lags),1);

    for i=1:length(lags)

        disp(i)
        smps_mean(:,i)=find_coincident_mean(bee_dataset.times-lags(i), smps_time, smps_sum, prof_len);


    end

    % save data
    lagged_smps=table();
    lagged_smps.times=bee_dataset.times;
    lagged_smps.bro_col=bee_dataset.bro_col;

    for i=1:length(lags)

        eval(['lagged_smps.smps_' num2str(i) 'h=smps_mean(:,i);']);

    end

    save lagged_smps.mat lagged_smps
    
end

%% get R^2

load('/home/kristof/work/BEEs/lagged_smps.mat')

bro_col=lagged_smps.bro_col;
lagged_smps=table2array(lagged_smps(:,3:end));

% % not good, since NaNs are in different places in each column -- discard a
% % lot of good data this way
% for i=2:72
%     lagged_smps(isnan(lagged_smps(:,i)),:)=[];
% end
% corr_coeffs=corr(lagged_smps);
    
R2_all=NaN(size(lagged_smps,2),1);
for i=1:size(lagged_smps,2)
    
    tmp2=lagged_smps(:,i);
    tmp1=bro_col(~isnan(tmp2));
    tmp2=tmp2(~isnan(tmp2));
    
    R2_all(i)=corr(tmp1,tmp2);
    
end

R2_all=R2_all.^2;

plot(R2_all)

xlabel('SMPS time lag (hours)')
ylabel('R^2 with BrO columns')

% end

