

%% load data
% for BrO:
    % vmr in ppm (prof, prof_err)
    % cnc in mu g /m^3
    % nd in moled/m^3 (prof_nd, prof_nd_err)

% for aerosol: extinction in km^-^1 (prof same as prof_nd)

% first value of each profile (alt=0) is the surface concentration

clf

prof_dir='/home/kristof/work/profile_retrievals/profile_results/eureka_';

year=2010;

tracegas=1;

% type='all'; % all profiles
type='filt'; % complete scans only

if tracegas
    
    % load data file
    load([prof_dir num2str(year) '/tracegas/profiles_' num2str(year) '_' type '.mat']);
    
    % convert to ppt (results are in ppm)
    prof=prof*1e6;
    prof_err=prof_err*1e6;
    
    prof(prof>30)=30;
    
else
    load([prof_dir num2str(year) '/aerosol/profiles_' num2str(year) '_' type '.mat']);
    prof(prof>10)=10;
end

[mean(info.DOFS), std(info.DOFS)]

%% plot profiles

% break up days and find gaps longer than 1h
prof_plot=prof;
ft_plot=ft;

ind_gap=find(ft(2:end)-ft(1:end-1)>1/24);

for i=1:length(ind_gap)

    % mean time of gap
    ft_plot=[ft_plot(1:ind_gap(i)); (ft_plot(ind_gap(i))+ft_plot(ind_gap(i)+1))/2;...
             ft_plot(ind_gap(i)+1:end)];
    
    % insert NaN profile
    prof_plot=[prof_plot(:,1:ind_gap(i)), NaN(21,1), prof_plot(:,ind_gap(i)+1:end)];
    
    % increment index
    ind_gap=ind_gap+1;
    
end

% plot color map

h=surf(ft_plot,alt,prof_plot,'EdgeColor','None', 'facecolor', 'interp'); hold on

view(2)

% colormap(jet(300))
colormap(hot(300))
c=colorbar;
grid off

xlim([ft_plot(1)-1; ft_plot(end)+1])

set(gca,'Color',[1 1 1]*0.6)

% part_prof=prof_nd(:,ind)*20000; % layers are 200m thick
% part_prof_err=prof_nd_err(:,ind)*20000;
% % errorbar(ft(ind)-subtract, sum(part_prof(1:top,:)),...
% %     (sum(part_prof_err(1:top,:).^2 ,1)).^0.5, 'b.','marker','square',...
% %     'markersize',4), hold on
% % errorbar(ft(ind)-subtract, sum(part_prof(top+1:end,:)),...
% %     (sum(part_prof_err(top+1:end,:).^2 ,1)).^0.5, 'r.','marker','square',...
% %     'markersize',4), hold on
% plot(ft(ind)-subtract, sum(part_prof(1:top,:)),'b.','marker','square',...
%     'markersize',4), hold on
% plot(ft(ind)-subtract, sum(part_prof(top+1:end,:)),'r.','marker','square',...
%     'markersize',4), hold on
% 
% legend(['0-' num2str(alt(top)+0.1) ' km'], [num2str(alt(top)+0.1) '-4 km'],...
%         'orientation','horizontal','location','northwest')



%%% delete all this below



% % select aerosol or tracegas
% option='a';
% if year==2016
%     %59 for leap years (2016), since data is already fractional day
%     subtract=59; 
% else
%     subtract=58;
% end
% 
% if option=='a'
%     filedir=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
%              num2str(year) '/aerosol/'];
% elseif option=='tg'
%     filedir=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
%              num2str(year) '/tracegas/'];
% end
% 
%    
% % profile file names
% search=[filedir '20*.mat' ];
% tmp = dir(search); 
% files={tmp.name}';
% 
% 
% if year==2017, limit=20; end % for aerosol extinction
% 
% tot_dof=[];
% 
% for i=1:size(files,2)
%     
%     % load file 
%     load([filedir,files{i}])
%     
%     if option=='a'
%         
%         % replace high values
%         if exist('limit','var')
%             ind=find(prof>limit);
%             prof(ind)=limit;
%         end
%         % only plot good profiles, filter by DoF
% %         ind=find(info(:,1)>1.5);
%         ind=find(info(:,1)>0);
%         
% %         % add zeros to each end to complete contours
% %         ft=[ft(1)-0.011;ft;ft(end)+0.011];
% %         fill=ones(size(prof(:,1)))*0.01;
% %         prof=[fill,prof,fill];
% %         ind=[1;ind+1;ind(end)+2];
% 
%     elseif option=='tg'
%         
%         % convert to ppt (results are in ppm)
%         prof=prof*1e6;
%         % only plot good profiles, filter by DoF
% %         ind=find(info(:,1)>0.8);
%         ind=find(info(:,1)>0);
% %         prof(prof>30)=30;
% 
%     end
% 
%         
%     %% plot profiles
%     figure(1)
% %     subplot(212)
% %     prof(prof>16)=16;
%     h=surf(ft(ind)-subtract,alt,prof(:,ind),'EdgeColor','None', 'facecolor', 'interp'); hold on
% 
%     figure(99)
%     h99_1=plot(ft(ind)-subtract,prof(4,ind),'ks'); hold on
% %     legend([h99_2,h99_1],{'Estimate','MAX-DOAS'},'Location','NorthWest')
%     
%     figure(50)
%     plot(ft(ind)-subtract,info(ind,1),'ks'), hold on
%     
%     % 0.5, 1.5, 3 for 2017
% 
% %     contour(ft(ind)-subtract,alt,prof(:,ind),'k:','levellist',[0.5]), hold on
% %     contour(ft(ind)-subtract,alt,prof(:,ind),'k--','levellist',[1.5]), hold on
% %     contour(ft(ind)-subtract,alt,prof(:,ind),'k-','levellist',[3.0]), hold on        
% 
% %     figure(2)
% %     errorbar(info(ind,4)-subtract, info(ind,2), info(ind,3), 'k','marker','square',...
% %         'linestyle','none','markersize',4), hold on
% %     plot(info(ind,4)-subtract, info(ind,2),'k.-','markersize',4), hold on
% 
%     
%     %% plot measured and modelled DSCDs
% %     figure(99)
% %     clist={'r.','g.','b.','y.','m.','c.','k.','b.'};
% % 
% %     nn=1;
% %     for elev=[-1,0,1,2,5,10,15,30]
% %         indelev=find(dscd(:,2)==elev);
% %         plot(ft_dscd(indelev),dscd(indelev,5),clist{nn},'markersize',7), hold on
% %         nn=nn+1;
% %     end
% %     
% %     nn=1;
% %     for elev=[-1,0,1,2,5,10,15,30]
% %         indelev=find(dscd(:,2)==elev);
% %         plot(ft_dscd(indelev),dscd(indelev,7),clist{nn},'marker','x'), hold on
% %         nn=nn+1;
% %     end
%     
%     %% plot VCD/extinction with/without errorbars
% %     figure(2)
% %     subplot(211)
% %     errorbar(info(ind,4)-subtract, info(ind,2), info(ind,3), 'k.-', 'markersize',8), hold on
%     
% %     figure(1)
% % %   
% % %     subplot(212)
% % %     plot(info(ind,4)-subtract, 100*info(ind,3)./info(ind,2), 'k.-', 'markersize',8), hold on
% % % 
%     if option=='tg'
% %         subplot(211)
%         top=4;
%         part_prof=prof_nd(:,ind)*20000; % layers are 200m thick
%         part_prof_err=prof_nd_err(:,ind)*20000;
% %         errorbar(ft(ind)-subtract, sum(part_prof(1:top,:)),...
% %             (sum(part_prof_err(1:top,:).^2 ,1)).^0.5, 'b.','marker','square',...
% %             'markersize',4), hold on
% %         errorbar(ft(ind)-subtract, sum(part_prof(top+1:end,:)),...
% %             (sum(part_prof_err(top+1:end,:).^2 ,1)).^0.5, 'r.','marker','square',...
% %             'markersize',4), hold on
%         plot(ft(ind)-subtract, sum(part_prof(1:top,:)),'b.','marker','square',...
%             'markersize',4), hold on
%         plot(ft(ind)-subtract, sum(part_prof(top+1:end,:)),'r.','marker','square',...
%             'markersize',4), hold on
%     
%         legend(['0-' num2str(alt(top)+0.1) ' km'], [num2str(alt(top)+0.1) '-4 km'],...
%                 'orientation','horizontal','location','northwest')
%         
% %         xlabel('Days of March, 2016 (UTC)')
%         ylabel('BrO VCD_{part} (molec/cm^2)')
% %         xlim([19.5, 23.49])
%         if year==2017, xlim([7, 21.5]), end
%         if year==2016, xlim([19, 23.5]), end
%         ylim([0,9.5]*1e13)
%         grid on
%         
%     end
% 
%     %% fit line to daily BrO VCD trend
% %     if option=='tg'
% %         % fit line to data
% %         [xData,yData] = prepareCurveData(info(ind,4)-subtract,info(ind,2));
% %         % fit line to data.
% %         [fitresult,~]=fit(xData,yData,'poly1');
% %         % plot result
% %         plot( info(ind,4)-subtract, fitresult(info(ind,4)-subtract), 'r--');
% %     end
%     
%     % plot DoFs
% %     figure(3)
% %     plot(info(ind,4)-subtract, info(ind,1), 'k.-'), hold on
%     
% %     tot_dof=[tot_dof; info(ind,1)];
% end
% 
% % figure(4)
% % histogram(tot_dof(1:67))
%     
% % figure(99)
% % legend('-1','0','1','2','5','10','15','30','modeled')
% 
% figure(1)
% % subplot(212)
% view(2)
% 
% ylim([0,4])
% xlim([7,21])
% 
% colormap(jet(300))
% c=colorbar;
% grid off
% 
% if option=='a'
%     ylabel(c,'Extinction (km^-^1)')
% elseif option=='tg'
% %     xlabel(c,'BrO VMR (pptv)')
%     xlabel(c,'MAX-DOAS BrO (pptv)')
% end
% 
% % xlabel(['Days of March, ' num2str(year) ' (UTC)'])
% ylabel('Altitude (km)')
% 
% % xlim([19.5,23.49])
% 
% % figure(2)
% % % subplot(212)
% % if option=='a'
% %     ylabel('Optical thickness')
% % elseif option=='tg'
% % %     ylabel('BrO VCD, 0-4 km (molec/cm^2)')
% %     ylabel('BrO VCD_{partial} (molec/cm^2)')
% % end
% % xlabel(['Days of March, ' num2str(year) ' (UTC)'])
% % % xlim([19.5,23.49])
% 
% % clearvars