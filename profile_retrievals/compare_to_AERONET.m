

%% load data

% load and merge aeronet data
aeronet_all=[];

load('/home/kristof/work/AERONET/PEARL_AOD_2015_5_2015_5.mat')
aeronet_all=[aeronet_all; aeronet];
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

% load retrieved profiles
load('/home/kristof/work/profile_retrievals/profile_results/aerosol_profiles_filt_all.mat')

ind_bad=find(info.DOFS<0.7 | info.col>5);

info(ind_bad,:)=[];
times(ind_bad,:)=[];

%% convert to AOD at 360 nm using aengstrom exponent
% HeiPro retrieval produces AOD at 360 nm
angstrom = - ( log(aeronet.AOD_340nm./aeronet.AOD_380nm) / log(340/380) );
aod_360=aeronet.AOD_380nm.*(360.8/380).^(-angstrom);

%% find coincident measurements, +-12 min
[ind_retr,ind_aeronet]=find_coincidences_time_general(times,aeronet.DateTime,0.2);
abs_diff=info.col(ind_retr)-aod_360(ind_aeronet);
rel_diff=(info.col(ind_retr)./aod_360(ind_aeronet)-1)*100;

figure
plot(times(ind_retr),info.col(ind_retr),'bo'), hold on
plot(aeronet.DateTime(ind_aeronet),aod_360(ind_aeronet),'k*')

figure
dscatter(aod_360(ind_aeronet),info.col(ind_retr))
xlim([0,1])
ylim([0,1])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

tmp=find(year(times(ind_retr))==2015);
ind_retr2=ind_retr(tmp(end)+1:end);
ind_aeronet2=ind_aeronet(tmp(end)+1:end);

figure
dscatter(aod_360(ind_aeronet2),info.col(ind_retr2))
xlim([0,1])
ylim([0,1])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

abs_diff2=info.col(ind_retr2)-aod_360(ind_aeronet2);
rel_diff2=(info.col(ind_retr2)./aod_360(ind_aeronet2)-1)*100;

disp([mean(abs_diff), mean(abs_diff2)])
disp([mean(rel_diff), mean(rel_diff2)])
