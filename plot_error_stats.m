
% CAMS error limits
%
% [UVVIS.DOAS.ZENITH.O3]
% O3.COLUMN.STRATOSPHERIC_SCATTER.SOLAR.ZENITH_UNCERTAINTY.SYSTEMATIC.STANDARD=[2.8,10]
% O3.COLUMN.STRATOSPHERIC_SCATTER.SOLAR.ZENITH_UNCERTAINTY.RANDOM.STANDARD=[3.5,5]
%
% [UVVIS.DOAS.ZENITH.NO2]
% NO2.COLUMN.STRATOSPHERIC_SCATTER.SOLAR.ZENITH_UNCERTAINTY.SYSTEMATIC.STANDARD=[1.5,24]
% NO2.COLUMN.STRATOSPHERIC_SCATTER.SOLAR.ZENITH_UNCERTAINTY.RANDOM.STANDARD=[2.8,18]



cd('/home/kristof/work/GBS/VCD_results')

msize=12;
xlim_max=360;
redist=3;

tg=1;
    
figure(1)
figure(2)

for year=[2014:2018]

    if tg==1
        load(['UT-GBS_O3_VCD_' num2str(year) '.mat'])
    else
        load(['UT-GBS_NO2_VCD_' num2str(year) '.mat'])
    end

    [ind,VCD_table]=filter_VCD_output(tg,VCD_table,rcd_S,1);

    VCD_table=VCD_table(ind,:);

    sys=(VCD_table.sigma_mean_vcd./VCD_table.mean_vcd)*100;
    rand=(VCD_table.std_vcd./VCD_table.mean_vcd)*100;
        
    if tg==1
        sys_rcd=sqrt(sys.^2 - 3.1^2);
    else
        sys_rcd=sqrt(sys.^2 - 9.43^2);
    end

    if redist==tg
        sys=sqrt(sys_rcd.^2 + 5^2);
        rand=sqrt(rand.^2 + 8^2);
    end
        
    figure(1)
%     subplot(311)
    subplot(321)
    plot(VCD_table.fd,sys,'.','markersize',msize), hold on

%     subplot(312)
%     plot(VCD_table.fd,rand,'.','markersize',msize), hold on

%     subplot(312)
    subplot(323)
    plot(VCD_table.fd,sys_rcd,'.','markersize',msize), hold on
    
%     subplot(313)
    subplot(325)
    plot(rcd_S.mean.day,rcd_S.mean.diff/2,'.','markersize',msize), hold on
    
    figure(2)
    subplot(221)
    plot(VCD_table.fd,rand,'.','markersize',msize), hold on

end

figure(1)

% subplot(311)
subplot(321)
grid on

ylabel('Systematic error (%)')
xlim([50,xlim_max])

if tg==1
    plot([50,xlim_max],[4,4],'k--','linewidth',1.3)
else
    plot([50,xlim_max],[9,9],'k--','linewidth',1.3)
    ylim([0,60])
end

legend('2014','2015','2016','2017','limit')

% subplot(312)
subplot(323)
grid on

% ylabel('Random error (%)')
ylabel('RCD error (%)')
xlim([50,xlim_max])

if tg==1
    plot([50,xlim_max],[4,4],'k--','linewidth',1.3)
else
    plot([50,xlim_max],[9,9],'k--','linewidth',1.3)
    ylim([0,60])
end

legend('2014','2015','2016','2017','limit')

% subplot(313)
subplot(325)
legend('2014','2015','2016','2017')
grid on

ylabel('\DeltaRCD/2 (molec/cm2)')
xlabel('Day of the year (UTC)')
xlim([50,xlim_max])

figure(2)

subplot(221)
grid on

ylabel('Random error (%)')
xlim([50,xlim_max])

if tg==1
    plot([50,xlim_max],[5,5],'k--','linewidth',1.3)
else
    plot([50,xlim_max],[12,12],'k--','linewidth',1.3)
%     ylim([0,60])
end

legend('2014','2015','2016','2017','limit')


%%%%%%%%%%%%%

if tg==1
    load('UT-GBS_O3_VCD_all.mat')
else
    load('UT-GBS_NO2_VCD_all.mat')
end

sys=(reanalysis.sigma_mean_vcd./reanalysis.mean_vcd)*100;
rand=(reanalysis.std_vcd./reanalysis.mean_vcd)*100;

if tg==1
    sys_rcd=sqrt(sys.^2 - 3.1^2);
else
    sys_rcd=sqrt(sys.^2 - 9.43^2);
end

if redist==tg
    sys=sqrt(sys_rcd.^2 + 5^2);
    rand=sqrt(rand.^2 + 8^2);
end

figure(1)

% subplot(211)
subplot(3,2,[2 4])
for i=60:20:280
    
    ind=find(reanalysis.day>=i & reanalysis.day<i+20);
%     ind=find(reanalysis.day>=i & reanalysis.day<i+20 & sys<=8);
    
    plot(i+10,mean(sys(ind)),'kx'), hold on
    plot(i+10,median(sys(ind)),'ko'), hold on

end

legend('mean','median')

ylabel('Systematic error for 20 day intervals (%)')
xlabel('Day of the year (UTC)')

grid on
grid minor


% subplot(212)
subplot(326)
if tg==1
    histogram(sys,50,'Normalization','probability','BinWidth',.25), hold on
    plot([median(sys),median(sys)],[0,.2],'r-','linewidth',1.3)
    plot([4,4],[0,.2],'k--','linewidth',1.3)
else
    histogram(sys,50,'Normalization','probability','BinWidth',1), hold on
    plot([median(sys),median(sys)],[0,.3],'r-','linewidth',1.3)
    plot([9,9],[0,.3],'k--','linewidth',1.3)
    if redist==tg
        xlim([0,50])
    else
        xlim([5,50])
    end
end

legend('data','median','QA/QC limit')

xlabel('Year-round systematic error (%)')
ylabel('Probability')




figure(2)

subplot(2,2,[2 4])

for i=60:20:280
    
    ind=find(reanalysis.day>=i & reanalysis.day<i+20);
%     ind=find(reanalysis.day>=i & reanalysis.day<i+20 & sys<=10);
    
    plot(i+10,mean(rand(ind)),'kx'), hold on
    plot(i+10,median(rand(ind)),'ko'), hold on

end

legend('mean','median')

ylabel('Random error for 20 day intervals (%)')
xlabel('Day of the year (UTC)')

grid on
grid minor


subplot(223)

if tg==1
    histogram(rand,50,'Normalization','probability','BinWidth',.25), hold on
    plot([median(rand),median(rand)],[0,.2],'r-','linewidth',1.3)
    plot([5,5],[0,.2],'k--','linewidth',1.3)
else
    histogram(rand,50,'Normalization','probability','BinWidth',1), hold on
    plot([median(rand),median(rand)],[0,.3],'r-','linewidth',1.3)
    plot([12,12],[0,.3],'k--','linewidth',1.3)
%     xlim([5,50])
end

legend('data','median','QA/QC limit')

xlabel('Year-round random error (%)')
ylabel('Probability')


