%% plot FYSI, MYSI, water, and land contact 
% as a function of BrO percentile and back trajectory length


c_list=[[1 0 0];...
        [0 1 0];...
        [0 0 1];...
        [0 1 1];...
        [1 0 1];...
        [.8 .8 0];...
        [.5 0 0];...
        [0 .5 0];...
        [0 0 .5];...
        [0 .5 .5]];

c_list=flipud(c_list);

%% gouped by wind direction
    
figure
set(gcf, 'Position', [100, 100, 1200, 650]);
fig_ax = tight_subplot(2,1,[0.07,0.01],[0.1,0.06],[0.02,0.16]);
% fig_ax = tight_subplot(2,1,[0.07,0.01],[0.08,0.06],[0.08,0.16]);

% select data to plot
plot_type_all={'FYSI','MYSI'};

for sp=1:2
    
    axes(fig_ax(sp))
    
    plot_type=plot_type_all{sp};

    switch plot_type
        case 'FYSI'
            box_data=bee_fp.fysi_3day;
        case 'MYSI'
            box_data=bee_fp.mysi_3day;
        case 'water'
            box_data=bee_fp.water_3day;
        case 'land'
            box_data=bee_fp.land_3day;
    end

    % create indices to plot bar plots with
    % each number in ascending order plots a different bar, code below sorts
    % BrO percentiles by wind direction, since wdir indices are 1-4  
    % (gives 100+ for N, 200+ for SE, 300+ for other)
    bee_fp.wdir(bee_fp.wdir==4)=3; % merge undetermined wdirs with 'other'
    box_group=bee_fp.bro_pc+100*bee_fp.wdir;
    % box_group=bee_fp.ssa_pc+100*bee_fp.wdir;

    % set positions to separate bar groups by wdir
    a=0.08; % controls width of each group
    nbars=9; % number of bars in each group (number of unique bro_pc indices)
    box_pos=sort([1:a:1+a*nbars,2:a:2+a*nbars,3:a:3+a*nbars]);

    boxplot(box_data,box_group,'positions',box_pos,'colors',c_list,'plotstyle','compact',...
            'symbol','s','jitter',0.3)
%     boxplot(box_data,box_group,'positions',box_pos,'colors',c_list,...
%             'symbol','s','jitter',0.3,'widths',0.06,'symbol','o')
%     set(findobj(gca,'type','line'),'linew',1.5)
        
    tick_pos=[];
    for i=1:3
        tick_pos=[tick_pos,mean(box_pos((nbars+1)*i-nbars:(nbars+1)*i))];
    end

    set(gca,'xtick',tick_pos)
    set(gca,'xticklabel',{'N','SE','other'})
    set(findall(gca,'tag','Outliers'),'MarkerSize',3);

    xlim([0.9,3.1+a*nbars])
    ylim([-0.7,10.7]*1e13)
    set(gca,'ygrid','on')

    set(findall(gca,'-property','FontSize'),'FontSize',17)

    if sp==1
        
        ll=legend(flipud(findall(gca,'Tag','Box')),...
               {'<10^t^h','10^t^h-20^t^h','20^t^h-30^t^h','30^t^h-40^t^h','40^t^h-50^t^h',...
                '50^t^h-60^t^h','60^t^h-70^t^h','70^t^h-80^t^h','80^t^h-90^t^h','>90^t^h'},...
               'Orientation','vertical','location','east');

        ll.FontSize=11;

        % [left bottom width height]
%         ll_pos=get(ll,'position');
        ll.Position(1)=0.865;
        ll.Position(2)=0.3;
        
        text(1.09,0.4,...
             sprintf('BrO column\npercentiles'), 'color','k','Units','normalized',...
             'fontsize',13,'HorizontalAlignment','center')
        
    else
        xlb=xlabel('Mean wind direction');
        xlb.Position(2)=xlb.Position(2)*1.9;
    end
    
    ylabel([plot_type ' contact (s m^2)'])
    
end



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


