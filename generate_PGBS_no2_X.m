%% generate x.xs for PGBS, for large NO2 window
load('/home/kristof/work/GBS/x_c-s_for_pgbs/resid_2007_d181-183');

% for day=183;
% % time window around noon in hours
% tw=2;
% 
% noon=0.7431; % 17:50 utc
% ind=find(fd>=day+noon-tw/24 & fd<=day+noon+tw/24);
% 
% figure
% plot(lambda,mean(resid(ind,:)),'kx'), hold on
% 
% tmp=boxcar(lambda,mean(resid(ind,:)),20);
% % plot(lambda,tmp,'linewidth',2)
% % plot(lambda,boxcar(lambda,mean(resid(ind,:)),40),'linewidth',2)
% 
% 
% x_xs=boxcar(lambda,tmp,10);
% plot(lambda,x_xs,'linewidth',2)
% 
% out=[lambda',x_xs'];
% dlmwrite('/home/kristof/work/GBS/x_c-s_for_pgbs/X_P1_resid_20_no2.xs',...
%          out,'delimiter','\t','precision',16)
% 
%      
% % X_P1_resid_20_no2.xs: d183 in 2007, 2h around noon, smoothed twice (20 and 10 pts)   
% end

%% compare 2007 UT-GBS and PEARL-GBS NO2 dSCDs
% like appendix B in Cristen's thesis, to see whether there is a
% correlation between NO2 dSCD and X

load('/home/kristof/work/GBS/x_c-s_for_pgbs/dSCDs_2007_ut+p');

no2diff=[];
no2x=[];
times=[];
dscds=[];
for i=1:617
    [tdiff,ind]=min(abs(data.Fractionalday - test.Fractionalday(i)));
    if tdiff<0.01 %14.4 min
        no2diff=[no2diff;test.NO2SlColno2(i) - data.NO2SlColno2(ind)];
        no2x=[no2x;test.NO2SlColX(i)];
        times=[times;test.Fractionalday(i)];
        dscds=[dscds;test.NO2SlColno2(i), data.NO2SlColno2(ind)];
    end
end

% plot(no2x,no2diff,'ko')
% plot(times,no2diff,'ko')

plot(dscds(:,1),dscds(:,2),'ko'), hold on
plot([-0.5,3.5]*1e16,[-0.5,3.5]*1e16,'k-')




