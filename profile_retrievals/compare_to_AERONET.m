

%% load data

% load and merge aeronet data
aeronet_all=[];

% load('/home/kristof/work/AERONET/PEARL_AOD_2015_5_2015_5.mat')
% aeronet_all=[aeronet_all; aeronet];
load('/home/kristof/work/AERONET/PEARL_AOD_2016_4_2016_5.mat')
aeronet_all=[aeronet_all; aeronet];
load('/home/kristof/work/AERONET/PEARL_AOD_2017_3_2017_5.mat')
aeronet_all=[aeronet_all; aeronet];
load('/home/kristof/work/AERONET/PEARL_AOD_2018_4_2018_5.mat')
aeronet_all=[aeronet_all; aeronet];

aeronet=aeronet_all;
clearvars aeronet_all

aeronet(aeronet.AOD_340nm==-999,:)=[];
aeronet(aeronet.AOD_380nm==-999,:)=[];

if 0
    % load retrieved profiles
    load('/home/kristof/work/profile_retrievals/profile_results/aerosol_profiles_filt_all.mat')

    % filter data
    ind_bad=find(times.Year==2015 & (info.DOFS<0.7 | info.col>5));

    info(ind_bad,:)=[];
    times(ind_bad,:)=[];
    
    time_data=times;
    ydata=info.col;
    
else
    
    % load/filter BEE dataset
    load('/home/kristof/work/BEEs/BEE_dataset_all.mat')

    % exclude 2015 (data already filtered)
    filt_ind=find(bee_dataset.times.Year==2015 & bee_dataset.aer_ext>1);

    bee_dataset(filt_ind,:)=[];
    
    time_data=bee_dataset.times;
    ydata=bee_dataset.aer_ext;

end
%% convert to AOD at 360 nm using aengstrom exponent
% HeiPro retrieval produces AOD at 360 nm
angstrom = - ( log(aeronet.AOD_340nm./aeronet.AOD_380nm) / log(340/380) );
aod_360=aeronet.AOD_380nm.*(360.8/380).^(-angstrom);

%% find coincident measurements, +-12 min
[ind_retr,ind_aeronet]=find_coincidences_time_general(time_data,aeronet.DateTime,0.2);
abs_diff=ydata(ind_retr)-aod_360(ind_aeronet);
rel_diff=(ydata(ind_retr)./aod_360(ind_aeronet)-1)*100;

figure
plot(time_data(ind_retr),ydata(ind_retr),'bo'), hold on
plot(aeronet.DateTime(ind_aeronet),aod_360(ind_aeronet),'k*')

figure
dscatter(aod_360(ind_aeronet),ydata(ind_retr))
xlim([0,1])
ylim([0,1])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% if filtering doesn't exclude 2015
% tmp=find(year(time_data(ind_retr))==2015);
% ind_retr2=ind_retr(tmp(end)+1:end);
% ind_aeronet2=ind_aeronet(tmp(end)+1:end);
% 
% figure
% dscatter(aod_360(ind_aeronet2),ydata(ind_retr2))
% xlim([0,1])
% ylim([0,1])
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% 
% abs_diff2=ydata(ind_retr2)-aod_360(ind_aeronet2);
% rel_diff2=(ydata(ind_retr2)./aod_360(ind_aeronet2)-1)*100;
% 
% disp([mean(abs_diff), mean(abs_diff2)])
% disp([mean(rel_diff), mean(rel_diff2)])
