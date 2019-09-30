% compare profiles (for the same day) with different settings

% need to be in correct aerosol or tracegas directory

% load both profiles
load 20160320_noadaptive_no0_1iter.mat
prof_1=prof; % reference profiles
vcd_1=info(:,2);
ft_1=info(:,4);

load 20160320_noadaptive_1iter_skip0.mat
prof_2=prof; % profiles to compare
vcd_2=info(:,2);
ft_2=info(:,4);

ft(end)=[];
ft_2(end)=[];
vcd_2(end)=[];
prof_2(:,end)=[];

% load 20160320_noadaptive-1.mat
% prof_2=prof(:,1:end-1); % profiles to compare
% vcd_2=info(1:end-1,2);
% ft_2=info(1:end-1,4);
% ft(end)=[];

clearvars -except alt ft prof_1 prof_2 vcd_1 vcd_2 ft_1 ft_2

% calculate differences
abs_diff=prof_2-prof_1;
rel_diff=abs_diff./prof_1;

abs_vcd=vcd_2-vcd_1;
rel_vcd=abs_vcd./vcd_1;

% plot results
figure(1)

%%% for general comparison
% surf(ft,alt,rel_diff,'EdgeColor','None', 'facecolor', 'interp')
% view(2)
% colormap(jet(300))
% colorbar
% xlim([ft(1),ft(end)])

%%% for comparing with and without -1 elevation
% plot results above (and including) instrument level
subplot(211)
surf(ft,alt(4:end-5),rel_diff(4:end-5,:),'EdgeColor','None', 'facecolor', 'interp')
view(2)
colormap(jet(300))
colorbar
xlim([ft(1),ft(end)])
ylim([0.6,3])

% plot results below instrument level
subplot(212)
surf(ft,alt(1:3),rel_diff(1:3,:),'EdgeColor','None', 'facecolor', 'interp')
view(2)
colormap(jet(300))
colorbar
xlim([ft(1),ft(end)])
ylim([0,0.4])

figure(2)
subplot(211)
plot(ft_1,vcd_1,'r.-'), hold on
plot(ft_2,vcd_2,'b.-'), hold on
subplot(212)
plot(ft_1, rel_vcd)
