%% plot FYSI, MYSI, water, and land contact 
% as a function of BrO percentile and back trajectory length

if 0
    load('/home/kristof/work/BEEs/BEE_dataset_flexpart.mat')
    for f=1:5 % back traj length in days

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/FP_FYSI_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.fysi_' num2str(f) 'day=FP_SI_contact;']);

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/FP_MYSI_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.mysi_' num2str(f) 'day=FP_SI_contact;']);

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/FP_water_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.water_' num2str(f) 'day=FP_SI_contact;']);

    end
end

c_list=[[1 0 0];...
        [0 1 0];...
        [0 0 1];...
        [0 1 1];...
        [1 0 1]];


%% gouped by wind direction
    
figure
box_data=bee_fp.fysi_3day;
% box_group=bee_fp.bro_pc+10*bee_fp.wdir;
box_group=bee_fp.ssa_pc+10*bee_fp.wdir;
box_group(box_group>40)=NaN;

a=0.14;
box_pos=sort([1:a:1+a*4,2:a:2+a*4,3:a:3+a*4]);

boxplot(box_data,box_group,'positions',box_pos,'colors',c_list,'plotstyle','compact',...
        'symbol','s','jitter',0.3)

tick_pos=[];
for i=1:3
    tick_pos=[tick_pos,mean(box_pos(5*i-4:5*i))];
end

set(gca,'xtick',tick_pos)
set(gca,'xticklabel',{'N','SE','other'})


set(findall(gca,'tag','Outliers'),'MarkerSize',2);

legend(flipud(findall(gca,'Tag','Box')), {'<25','25-50','50-75','75-90','>90'},...
       'Orientation','horizontal','location','northwest');

set(gca,'ygrid','on')

xlabel('Mean wind direction')
ylabel('FYSI contact (s m^2)')

%% grouped by back traj length, 5 BrO percentile ranges for each

% figure
% % box_data=[bee_fp.mysi_2day;bee_fp.mysi_3day;bee_fp.mysi_4day];
% box_data=[bee_fp.fysi_2day;bee_fp.fysi_3day;bee_fp.fysi_4day];
% box_group=[bee_fp.bro_pc;bee_fp.bro_pc+10;bee_fp.bro_pc+20];
% 
% a=0.14;
% box_pos=sort([1:a:1+a*4,2:a:2+a*4,3:a:3+a*4]);
% 
% boxplot(box_data,box_group,'positions',box_pos,'colors',c_list,'plotstyle','compact',...
%         'symbol','s','jitter',0.3)
% 
% tick_pos=[];
% for i=1:3
%     tick_pos=[tick_pos,mean(box_pos(5*i-4:5*i))];
% end
% 
% set(gca,'xtick',tick_pos)
% set(gca,'xticklabel',{'2 days','3 days','4 days'})
% 
% 
% set(findall(gca,'tag','Outliers'),'MarkerSize',2);
% 
% legend(flipud(findall(gca,'Tag','Box')), {'<25','25-50','50-75','75-90','>90'},...
%        'Orientation','horizontal','location','northwest');
% 
% set(gca,'ygrid','on')
% 
% xlabel('Back trajectory duration')
% ylabel('FYSI contact (s m^2)')

%% grouped by BrO, 3 back traj lengths for each

% figure
% box_data=[bee_fp.fysi_2day;bee_fp.fysi_3day;bee_fp.fysi_4day];
% tmp=[bee_fp.bro_pc;bee_fp.bro_pc+10;bee_fp.bro_pc+20];
% 
% box_pos=sort([1:5,1.2:5.2,1.4:5.4]);
% 
% num=1;
% box_group=NaN(size(tmp));
% for i=1:5
% 
%     box_group(tmp==i)=num;
%     box_group(tmp==i+10)=num+1;
%     box_group(tmp==i+20)=num+2;
%     num=num+3;
%     
% end
% 
% boxplot(box_data,box_group,'positions',box_pos,'colors','rgb','plotstyle','compact')
% 
% tick_pos=[];
% for i=1:5
%     tick_pos=[tick_pos,mean(box_pos(3*i-2:3*i))];
% end
% 
% set(gca,'xtick',tick_pos)
% set(gca,'xticklabel',{'<25','25-50','50-75','75-90','>90'})
% 
% hLegend = legend(findall(gca,'Tag','Box'), {'4 day','3 day','2 day'});grid on

%% individual plots

% figure
% subplot(221)
% boxplot(bee_fp.fysi_3day,bee_fp.bro_pc,'labels',{'<25','25-50','50-75','75-90','>90'})
% 
% xlabel('BrO column percentile')
% ylabel('3 day FYSI contact (s m^2)')
% 
% % figure
% subplot(223)
% boxplot(bee_fp.mysi_3day,bee_fp.bro_pc,'labels',{'<25','25-50','50-75','75-90','>90'})
% 
% xlabel('BrO column percentile')
% ylabel('3 day MYSI contact (s m^2)')
% 
% % figure
% subplot(222)
% boxplot(bee_fp.fysi_5day,bee_fp.bro_pc,'labels',{'<25','25-50','50-75','75-90','>90'})
% 
% xlabel('BrO column percentile')
% ylabel('5 day FYSI contact (s m^2)')
% 
% % figure
% subplot(224)
% boxplot(bee_fp.mysi_5day,bee_fp.bro_pc,'labels',{'<25','25-50','50-75','75-90','>90'})
% 
% xlabel('BrO column percentile')
% ylabel('5 day MYSI contact (s m^2)')


