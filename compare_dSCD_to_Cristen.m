

% Cristen's results
load('/home/kristof/work/GBS/VCD_results/June2011_Reanalysis_All.mat')
% load other .mat file manually

%% compare dSCD values those found in Cristen's deluge folder
% % % 
% % % indc=find(cristen(:,7)<91 & cristen(:,12)<0.002 & cristen(:,21)>0);
% % % indme=find(qdoas_filt(:,7)<91 & qdoas_filt(:,181)>0);
% % % 
% % % % plot(cristen(indc,3),cristen(indc,21),'ro'); hold on
% % % % plot(qdoas_filt(indme,3),qdoas_filt(indme,181),'bx');
% % % 
% % % % find coincidences
% % % 
% % % cristen2=cristen(indc,:);
% % % me=qdoas_filt(indme,:);
% % % 
% % % [~,ime,ic]=intersect(me(:,3),cristen2(:,3));
% % % 
% % % 
% % % rel_diff=((me(ime,181)./cristen2(ic,21))-1)*100;
% % % time=me(ime,3);
% % % 
% % % 
% % % % mean(rel_diff)
% % % % std(rel_diff)
% % % 
% % % figure
% % % plot(time,rel_diff,'ko'), hold on


%% Compare langley plots
% % % 
% % % day=257;
% % % ampm=1;
% % % 
% % % 
% % % ind1=find(O3_t1_rcd.P2007.day==day & O3_t1_rcd.P2007.ampm==ampm);
% % % ind2=find(rcd_S.day==day & rcd_S.ampm==ampm);
% % % 
% % % % amf=[1.9,4.2]; % summer
% % % amf=[10,20]; % spring/fall
% % % 
% % % figure
% % % plot(amf,O3_t1_rcd.P2007.y_int(ind1)+O3_t1_rcd.P2007.m(ind1)*amf,'k-'); hold on
% % % plot(amf,rcd_S.y_int(ind2)+rcd_S.m(ind2)*amf,'r--'); hold on
% % % xlabel('Arbitrary AMF range')
% % % ylabel('DSCD from Langley fit parameters')

%% Compare RCDs
% figure
% plot(O3_t1_rcd.P2007.mean.day,O3_t1_rcd.P2007.mean.rcd,'ko','linewidth',1.5), hold on
% plot(rcd_S.mean.day,rcd_S.mean.rcd,'kx','linewidth',1.5), hold on
% xlabel('Day of the year')
% ylabel('RCD')

