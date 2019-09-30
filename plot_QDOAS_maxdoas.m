% to plot maxdoas results for PEARL-GBS

% load('/home/kristof/work/GBS/PEARL-GBS/2016/MAX-DOAS/maxdoas_bro.mat')
% load('/home/kristof/work/GBS/PEARL-GBS/2017/MAX-DOAS/maxdoas_bro.mat')

% figure('Position', [100, 100, 1410, 1050])
figure(1)

subtract=60; % 60 for leap years since data is in day of year
% subtract=59; % 59 for non-leap years since data is in day of year

% subtract=-1; % for date on x axis -- datetick converts DOY, but assumes leap year

day1 = 65; 
day2 = 131; 
% day1 = 66;
% day2 = 70;

indm1=find(dscd_S_m1.day >= day1 & dscd_S_m1.day <= day2);
ind0=find(dscd_S_0.day >= day1 & dscd_S_0.day <= day2);
ind1=find(dscd_S_1.day >= day1 & dscd_S_1.day <= day2);
ind2=find(dscd_S_2.day >= day1 & dscd_S_2.day <= day2);
ind5=find(dscd_S_5.day >= day1 & dscd_S_5.day <= day2);
ind10=find(dscd_S_10.day >= day1 & dscd_S_10.day <= day2);
ind15=find(dscd_S_15.day >= day1 & dscd_S_15.day <= day2);
ind30=find(dscd_S_30.day >= day1 & dscd_S_30.day <= day2);

msize=16;

% ax=subplot(2,1,1);
box on
hold on;
grid on
% grid minor
% ax.GridAlpha=0.4;
% grid minor
plot(dscd_S_m1.fd(indm1)-subtract,dscd_S_m1.mol_dscd(indm1),'r.','markersize', msize);
plot(dscd_S_0.fd(ind0)-subtract,dscd_S_0.mol_dscd(ind0),'g.','markersize', msize);
plot(dscd_S_1.fd(ind1)-subtract,dscd_S_1.mol_dscd(ind1),'b.','markersize', msize);
plot(dscd_S_2.fd(ind2)-subtract,dscd_S_2.mol_dscd(ind2),'y.','markersize', msize);
plot(dscd_S_5.fd(ind5)-subtract,dscd_S_5.mol_dscd(ind5),'m.','markersize', msize);
plot(dscd_S_10.fd(ind10)-subtract,dscd_S_10.mol_dscd(ind10),'c.','markersize', msize);
plot(dscd_S_15.fd(ind15)-subtract,dscd_S_15.mol_dscd(ind15),'k.','markersize', msize);
plot(dscd_S_30.fd(ind30)-subtract,dscd_S_30.mol_dscd(ind30),'ko','markersize', 3);
ylim([-1e14 10e14]);
% xlim([day1-subtract,day2+1.5-subtract])
xlim([day1-subtract,day2-subtract])
% xlabel('Fractional day')
ylabel('BrO DSCD (mol/cm^2)')
legend('-1\circ','0\circ','1\circ','2\circ','5\circ',...
    '10\circ','15\circ','30\circ','location','northeast','Orientation','horizontal')

% datetick('x','mmmdd','keeplimits')

% subplot(3,1,2);
% hold on;
% grid on
% grid minor
% plot(dscd_S_m1.fd(indm1),dscd_S_m1.rms(indm1),'r.');
% plot(dscd_S_0.fd(ind0),dscd_S_0.rms(ind0),'g.');
% plot(dscd_S_1.fd(ind1),dscd_S_1.rms(ind1),'b.');
% plot(dscd_S_2.fd(ind2),dscd_S_2.rms(ind2),'y.');
% plot(dscd_S_5.fd(ind5),dscd_S_5.rms(ind5),'m.');
% plot(dscd_S_10.fd(ind10),dscd_S_10.rms(ind10),'c.');
% plot(dscd_S_15.fd(ind15),dscd_S_15.rms(ind15),'k.');
% plot(dscd_S_30.fd(ind30),dscd_S_30.rms(ind30),'ko','markersize', 2);
% xlabel('Fractional day')
% ylabel('BrO RMS')
% ylim([0,0.0035])
% xlim([day1,day2+1.5])
% legend('-1','0','1','2','5','10','15','30','location','northeastoutside')

% % ax=subplot(2,1,2);
% % hold on;
% % grid on
% % % ax.GridAlpha=0.4;
% % % grid minor
% % plot(dscd_S_m1.fd(indm1)-subtract,dscd_S_m1.mol_dscd_o4(indm1)*1e40,'r.','markersize', msize);
% % plot(dscd_S_0.fd(ind0)-subtract,dscd_S_0.mol_dscd_o4(ind0)*1e40,'g.','markersize', msize);
% % plot(dscd_S_1.fd(ind1)-subtract,dscd_S_1.mol_dscd_o4(ind1)*1e40,'b.','markersize', msize);
% % plot(dscd_S_2.fd(ind2)-subtract,dscd_S_2.mol_dscd_o4(ind2)*1e40,'y.','markersize', msize);
% % plot(dscd_S_5.fd(ind5)-subtract,dscd_S_5.mol_dscd_o4(ind5)*1e40,'m.','markersize', msize);
% % plot(dscd_S_10.fd(ind10)-subtract,dscd_S_10.mol_dscd_o4(ind10)*1e40,'c.','markersize', msize);
% % plot(dscd_S_15.fd(ind15)-subtract,dscd_S_15.mol_dscd_o4(ind15)*1e40,'k.','markersize', msize);
% % plot(dscd_S_30.fd(ind30)-subtract,dscd_S_30.mol_dscd_o4(ind30)*1e40,'ko','markersize', 3);
% % ylim([-1000,7000]*1e40);
% % xlabel('Days of March, 2016 (UTC)')
% % ylabel('O_4 DSCD (mol^2/cm^5)')
% % xlim([day1-subtract,day2+1.5-subtract])
% % % legend('-1','0','1','2','5','10','15','30','location','northeastoutside')

% set font on plots
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman') 
% 
% f=gcf; 
% figpos=getpixelposition(f); 
% resolution=get(0,'ScreenPixelsPerInch'); 
% set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
% path='/home/kristof/work/summer_school/poster/'; 
% name='BrO_march_19-22'; 
% print(f,fullfile(path,name),'-dpng','-r300','-opengl') %save file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% For CMOS poster: plot march-april data for 2015-2017
% % 
% % years=[2015,2016,2017];
% % 
% % figure(1)
% % 
% % for i=1:3
% % 
% %     load(['/home/kristof/work/GBS/PEARL-GBS/' num2str(years(i)) '/MAX-DOAS/maxdoas_bro.mat'])
% %     
% %     if years(i)==2015
% %         dscd_S_2.day=80;
% %         dscd_S_2.fd=80.6;
% %         dscd_S_2.mol_dscd=0;
% %     end
% % 
% %     day1 = 60; 
% %     day2 = 120; 
% % 
% %     fdm1_avg=NaN(1,day2-day2+1);
% %     fd0_avg=NaN(1,day2-day2+1);
% %     fd1_avg=NaN(1,day2-day2+1);
% %     fd5_avg=NaN(1,day2-day2+1);
% %     
% %     dscd_m1_avg=NaN(1,day2-day2+1);
% %     dscd_0_avg=NaN(1,day2-day2+1);
% %     dscd_1_avg=NaN(1,day2-day2+1);
% %     dscd_5_avg=NaN(1,day2-day2+1);
% %     
% %     count=1;
% %     for j=[day1:day2]
% %         
% %         if j==98 && i==1, continue, end
% %         
% %         indm1=find(dscd_S_m1.day == j);
% %         ind0=find(dscd_S_0.day == j);
% %         ind1=find(dscd_S_1.day == j);
% %         ind5=find(dscd_S_5.day == j);
% %         
% %         fdm1_avg(count)=mean(dscd_S_m1.fd(indm1));
% %         dscd_m1_avg(count)=mean(dscd_S_m1.mol_dscd(indm1));
% % 
% %         fd0_avg(count)=mean(dscd_S_0.fd(ind0));
% %         dscd_0_avg(count)=mean(dscd_S_0.mol_dscd(ind0));
% % 
% %         fd1_avg(count)=mean(dscd_S_1.fd(ind1));
% %         dscd_1_avg(count)=mean(dscd_S_1.mol_dscd(ind1));
% % 
% %         fd5_avg(count)=mean(dscd_S_5.fd(ind5));
% %         dscd_5_avg(count)=mean(dscd_S_5.mol_dscd(ind5));
% % 
% %         count=count+1;
% %         
% %     end
% % 
% %     msize=24;
% % 
% %     ax=subplot(3,1,i);
% %     box on
% %     hold on;
% %     grid on
% %     % grid minor
% %     % ax.GridAlpha=0.4;
% %     % grid minor
% %     plot(fd5_avg,dscd_5_avg,'c.','markersize', msize);
% %     plot(fd1_avg,dscd_1_avg,'b.','markersize', msize);
% %     plot(fd0_avg,dscd_0_avg,'g.','markersize', msize);
% %     plot(fdm1_avg,dscd_m1_avg,'r.','markersize', msize);
% % 
% %     ylim([-1e14 8e14]);
% %     xlim([day1,day2+1])
% %     
% %     legend('5\circ','1\circ','0\circ','-1\circ','location','northeast','Orientation','horizontal')
% % 
% %     text(83,5e14,num2str(years(i)))
% %     
% %     if i==2
% %         hh=ylabel('Average BrO DSCD (molec/cm^2)');
% %         set(hh,'position',[57,3.5e14,-1])
% %     end
% %     if i==3, xlabel('Fractional day (UTC)'), end
% % 
% % end
