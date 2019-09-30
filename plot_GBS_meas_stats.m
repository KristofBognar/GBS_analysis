% plot number of measurements for spring campaign


%% load raw QDOAS files

year=2019;
% instr='UT-GBS';
instr='PEARL-GBS';

cur_dir=pwd();
qdoas_dir='/home/kristof/work/GBS/QDOAS_results/NDACC_RD_tables/';

% make list of QDOAS files
temp = dir([qdoas_dir instr '_' num2str(year) '*.mat']); 
f_list = {temp.name}; % cell array of file names

% loop over QDOAS files
cd(qdoas_dir)

data_all=[];


for file=f_list
    
    load(file{1})
    
    data_all=[data_all;data];

end

bin_edges=floor(data_all.Fractionalday(1))-0.5:floor(data_all.Fractionalday(end))+0.5;

%% ZS measurements
figure(1)
subplot(211)

% total number
histogram(floor(data_all.Fractionalday),bin_edges), hold on

% filter by SZA range
ind=find(data_all.SZA<=91 & data_all.SZA>86);
histogram(floor(data_all.Fractionalday(ind)), bin_edges)

xlabel(['Day of the year, ' num2str(year) ' (UTC)'])
ylabel('N.o. measurements')
legend('SZA=all','SZA=[86,91]','location','best')
title('Number of Zenith-sky measurements')
xlim([43,121])

%% MAX-DOAS measurements

% dunno how I counted last year
% count by difference in meas. fd (anything larger than 0.005 likely is a
% maxdoas scan in between) 
subplot(212)

diffs=data_all.Fractionalday(1:end-1)-data_all.Fractionalday(2:end);
ind=find(abs(diffs)>0.005);
histogram(floor(data_all.Fractionalday(ind)), bin_edges);
ylim([3,130])
xlim([43,121])

xlabel(['Day of the year, ' num2str(year) ' (UTC)'])
ylabel('N.o. scans')

if strcmp(instr,'UT-GBS')
    title('Number of MAX-DOAS scans (5 elevations each)')
else
    title('Number of MAX-DOAS scans (8 elevations each)')
end

%% set ffigure and font size
set(gcf, 'Position', [100, 100, 900, 550]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

f=gcf; 
f.Units = 'pixels';

figpos=getpixelposition(f); 
resolution=get(0,'ScreenPixelsPerInch'); 
set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 

% path=['/home/kristof/work/campaigns/campaign_' num2str(year) '/post_campaign']; 
% name=[instr '_stats'];
% print(f,fullfile(path,name),'-dpng','-r300','-opengl') %save file


cd(cur_dir)


% 2017 version

% % %% UT-GBS
% % %% load raw QDOAS files
% % load /home/kristof/work/GBS/UT-GBS/2017/VCD/ozone/ozone_v2_86-90.mat
% % 
% % bin_edges=dscd_S.day(1)-0.5:dscd_S.day(end)+0.5;
% % 
% % %% ZS measurements
% % figure(1)
% % subplot(211)
% % 
% % % total number
% % histogram(dscd_S.day,bin_edges), hold on
% % 
% % % filter by SZA range
% % ind=find(dscd_S.sza<=91 & dscd_S.sza>86);
% % histogram(dscd_S.day(ind), bin_edges)
% % 
% % xlabel('Day of the year, 2017 (UTC)')
% % ylabel('N.o. measurements')
% % legend('SZA=all','SZA=[86,91]','location','best')
% % title('Number of Zenith-sky measurements')
% % xlim([45,115])
% % 
% % %% MAX-DOAS measurements
% % 
% % % dunno how I counted last year
% % % count by difference in meas. fd (anything larger than 0.005 likely is a
% % % maxdoas scan in between) 
% % subplot(212)
% % 
% % diffs=dscd_S.fd(1:end-1)-dscd_S.fd(2:end);
% % ind=find(abs(diffs)>0.005);
% % histogram(dscd_S.day(ind), bin_edges);
% % ylim([5,100])
% % xlim([45,115])
% % 
% % xlabel('Day of the year, 2017 (UTC)')
% % ylabel('N.o. scans')
% % title('Number of MAX-DOAS scans (5 elevations each)')
% % 
% % %% PEARL-GBS
% % %% load raw QDOAS files
% % load /home/kristof/work/GBS/PEARL-GBS/2017/VCD/no2/no2_v1_86_91.mat
% % 
% % bin_edges=dscd_S.day(1)-0.5:dscd_S.day(end)+0.5;
% % 
% % %% ZS measurements
% % figure(2)
% % subplot(211)
% % 
% % % total number
% % histogram(dscd_S.day,bin_edges), hold on
% % 
% % % filter by SZA range
% % ind=find(dscd_S.sza<=91 & dscd_S.sza>86);
% % histogram(dscd_S.day(ind), bin_edges)
% % 
% % ylabel('N.o. measurements')
% % legend('SZA=all','SZA=[86,91]','location','best')
% % title('Number of Zenith-sky measurements')
% % xlim([59,115])
% % ylim([0,1200])
% % 
% % %% MAX-DOAS measurements
% % 
% % load /home/kristof/work/GBS/PEARL-GBS/2017/MAX-DOAS/maxdoas_bro.mat
% % 
% % subplot(212)
% % histogram(dscd_S_30.day,bin_edges)
% % 
% % ylabel('N.o. scans')
% % xlabel('Day of the year, 2017 (UTC)')
% % title('Number of MAX-DOAS scans (8 elevations each)')
% % xlim([59,115])
