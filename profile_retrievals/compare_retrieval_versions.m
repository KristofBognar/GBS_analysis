% compare Eureka aerosol retrievals with different number of iteration steps
%
% 5 iteration steps works best: it's closer to dSCDs (with slightly more
% scatter), has slightly higher DOFS, and compares much better to AERONET
% than 10 or 10 iter versions. 10/15 iterations also lead to unphysical
% extinction peaks above 3km.

%% load data versions

load('/home/kristof/work/AERONET/PEARL_AOD_2018_4_2018_5.mat')
angstrom = - ( log(aeronet.AOD_340nm./aeronet.AOD_380nm) / log(340/380) );
aod_360=aeronet.AOD_380nm.*(360.8/380).^(-angstrom);

% test version 1 (10 iter steps)
load('/home/kristof/work/profile_retrievals/profile_results/eureka_2018_aer_10iter/aerosol/profiles_2018_filt.mat');
dscd(dscd.elev==90,:)=[];
dscd(dscd.retr==0,:)=[];

times_10=times;
dscd_10=dscd;
info_10=info;
prof_10=prof;

% test version 2 (15 iter steps)
load('/home/kristof/work/profile_retrievals/profile_results/eureka_2018_aer_15iter/aerosol/profiles_2018_filt.mat');
dscd(dscd.elev==90,:)=[];
dscd(dscd.retr==0,:)=[];

times_15=times;
dscd_15=dscd;
info_15=info;
prof_15=prof;

% original retrieval (5 iter steps)
load('/home/kristof/work/profile_retrievals/profile_results/eureka_2018/aerosol/profiles_2018_filt.mat');
dscd(dscd.elev==90,:)=[];
dscd(dscd.O4retr==0,:)=[];

times_orig=times;
dscd_orig=dscd;
info_orig=info;
prof_orig=prof;

dscd_orig(dscd_orig.date_time<dscd_10.date_time(1) | ...
          dscd_orig.date_time>dscd_10.date_time(end),:)=[];
ind=find(times_orig<times_10(1) | times_orig>times_10(end));
times_orig(ind)=[];
info_orig(ind,:)=[];
prof_orig(:,ind)=[];

%% compare

% mean(dscd_orig.O4retr-dscd_orig.O4meas)
% mean(dscd_10.retr-dscd_10.meas)
% mean(dscd_15.retr-dscd_15.meas)
% std(dscd_orig.O4retr-dscd_orig.O4meas)
% std(dscd_10.retr-dscd_10.meas)
% std(dscd_15.retr-dscd_15.meas)

% mean((dscd_orig.O4retr./dscd_orig.O4meas)-1)*100
% mean((dscd_10.retr./dscd_10.meas)-1)*100
% mean((dscd_15.retr./dscd_15.meas)-1)*100
% std((dscd_orig.O4retr./dscd_orig.O4meas)-1)*100
% std((dscd_10.retr./dscd_10.meas)-1)*100
% std((dscd_15.retr./dscd_15.meas)-1)*100

% mean(info_orig.DOFS)
% mean(info_10.DOFS)
% mean(info_15.DOFS)

% figure, hold on
% plot(aeronet.DateTime,aod_360,'go')
% plot(times_orig,info_orig.col,'ko')
% plot(times_10,info_10.col,'rx')
% plot(times_15,info_15.col,'b+')


% figure, hold on
% plot(dscd_orig.date_time,dscd_orig.O4retr-dscd_orig.O4meas,'ko')
% plot(dscd_10.date_time,dscd_10.retr-dscd_10.meas,'rx')
% plot(dscd_15.date_time,dscd_15.retr-dscd_15.meas,'b+')


% [ind_orig,aeronet_orig]=find_coincidences_time_general(times_orig,aeronet.DateTime,0.25);
% abs_diff_orig=info_orig.col(ind_orig)-aod_360(aeronet_orig);
% rel_diff_orig=(info_orig.col(ind_orig)./aod_360(aeronet_orig)-1)*100;
% 
% [ind_10,aeronet_10]=find_coincidences_time_general(times_10,aeronet.DateTime,0.25);
% abs_diff_10=info_10.col(ind_10)-aod_360(aeronet_10);
% rel_diff_10=(info_10.col(ind_10)./aod_360(aeronet_10)-1)*100;
% 
% [ind_15,aeronet_15]=find_coincidences_time_general(times_15,aeronet.DateTime,0.25);
% abs_diff_15=info_15.col(ind_15)-aod_360(aeronet_15);
% rel_diff_15=(info_15.col(ind_15)./aod_360(aeronet_15)-1)*100;

% disp([mean(abs_diff_orig),std(abs_diff_orig)])
% disp([mean(abs_diff_10),std(abs_diff_10)])
% disp([mean(abs_diff_15),std(abs_diff_15)])
% 
% disp([mean(rel_diff_orig),std(rel_diff_orig)])
% disp([mean(rel_diff_10),std(rel_diff_10)])
% disp([mean(rel_diff_15),std(rel_diff_15)])

% [~,~,r1]=line_fit(aod_360(aeronet_orig),info_orig.col(ind_orig))
% [~,~,r2]=line_fit(aod_360(aeronet_10),info_10.col(ind_10))
% [~,~,r3]=line_fit(aod_360(aeronet_15),info_15.col(ind_15))

% plot(aod_360(aeronet_orig),info_orig.col(ind_orig),'k.'), hold on
% plot(aod_360(aeronet_10),info_10.col(ind_10),'r.')
% plot(aod_360(aeronet_15),info_15.col(ind_15),'b.')


%%%%%%%%%%%%%%%%%%%%%%%

prof_plot=prof_orig;
ft_plot=fracdate(times_orig);

prof_plot(prof_plot>5)=5;
prof_plot(end,1)=3e-4;

% ind_gap=find(ft_plot(2:end)-ft_plot(1:end-1)>1/24);
% for i=1:length(ind_gap)
%     % mean time of gap
%     ft_plot=[ft_plot(1:ind_gap(i)); (ft_plot(ind_gap(i))+ft_plot(ind_gap(i)+1))/2;...
%              ft_plot(ind_gap(i)+1:end)];
%     % insert NaN profile
%     prof_plot=[prof_plot(:,1:ind_gap(i)), NaN(21,1), prof_plot(:,ind_gap(i)+1:end)];
%     % increment index
%     ind_gap=ind_gap+1;
% end

% % plot color map
% figure
% h=surf(ft_plot,alt,log10(prof_plot),'EdgeColor','None', 'facecolor', 'interp'); hold on
% colorbar();
% view(2)
% 
% %%%
% figure
% prof_plot=prof_10;
% ft_plot=fracdate(times_10);
% 
% prof_plot(prof_plot>5)=5;
% prof_plot(end,1)=3e-4;
% 
% h=surf(ft_plot,alt,log10(prof_plot),'EdgeColor','None', 'facecolor', 'interp'); hold on
% colorbar();
% view(2)

% %%%
% figure
% prof_plot=prof_15;
% ft_plot=fracdate(times_15);
% 
% prof_plot(prof_plot>10)=10;
% prof_plot(end,1)=3e-4;
% 
% h=surf(ft_plot,alt,prof_plot,'EdgeColor','None', 'facecolor', 'interp'); hold on
% colorbar();
% view(2)

[ind1,ind2]=find_coincidences_time_general(times_10,times_orig,0.25);

plot(times_10(ind1),((info_10.col(ind1)./info_orig.col(ind2))-1)*100,'ko')
