%%% plot yearly ozone and NO2 in the spring, to look at ozone depletion
%%% plot modelled on Adams et al., 2012

du=2.687e16;

curr_yr=2020;

o3_only=1;

instr='UT-GBS';
% instr='PEARL-GBS';

highlight=[2000,2005,2007,2011,2020];
x_lim=[50,112];

plot_colors=flipud([[1 0 0];...
             [0 0 1];...
             [0 1 0];...
             [1 1 0];...
             [1 0 1];...
                    ]);

% plot_colors=flipud(jet(length(highlight)));

% merged GBS data, up to 2017 (from merged VCDs, already filtered)
load('/home/kristof/work/satellite_validation/GBS_NO2.mat')
load('/home/kristof/work/satellite_validation/GBS_O3.mat')

%% add later ozone data
% load merged files only, since those are filtered
load('/home/kristof/work/GBS/VCD_results/UT-GBS_O3_VCD_all.mat')
gbs_o3=[gbs_o3; reanalysis(reanalysis.year>2017,:)];

% check if some data is still RD (no standard processing yet)
check=max(reanalysis.year)+1:curr_yr;

for yr=check
    
    load(['/home/kristof/work/GBS/VCD_results/NDACC_RD/UT-GBS_O3_VCD_' num2str(yr) 'all.mat'])
    gbs_o3=[gbs_o3; reanalysis];
    
end

% convert to DU
gbs_o3.mean_vcd=gbs_o3.mean_vcd./du;

% remove extra data
gbs_o3(gbs_o3.fractional_time>x_lim(2),:)=[];

%% add later NO2 data
% load merged files only, since those are filtered
load('/home/kristof/work/GBS/VCD_results/UT-GBS_NO2_VCD_all.mat')
gbs_no2=[gbs_no2; reanalysis(reanalysis.year>2017,:)];

% check if some data is still RD (no standard processing yet)
check=max(reanalysis.year)+1:curr_yr;

for yr=check
    
    load(['/home/kristof/work/GBS/VCD_results/NDACC_RD/UT-GBS_NO2_VCD_' num2str(yr) 'all.mat'])
    gbs_no2=[gbs_no2; reanalysis];
    
end

% remove extra data
gbs_no2(gbs_no2.fractional_time>x_lim(2),:)=[];

%% plot results

if ~o3_only

    figure
    set(gcf, 'Position', [100, 100, 1000, 600]);
    fig_ax = tight_subplot(2,1,[0.07,0.07],[0.1,0.08],[0.08,0.13]);

    txt_y=0.97;

    for yr=unique(gbs_o3.year)'


        ind=find(highlight==yr);

        if ~isempty(ind)
            plot_c=plot_colors(ind,:);
        else
            plot_c=[.5 .5 .5];
        end

        %%% oxone plot
        axes(fig_ax(1))

        ind=gbs_o3.year==yr;
        xx=mjd2k_to_date(gbs_o3.mjd2k);
        xx.Year=0;

        plot(xx(ind),gbs_o3.mean_vcd(ind),'s','color',plot_c,...
             'linewidth',1.2), hold on

        text(1.05,txt_y,num2str(yr),'units','normalized','fontweight','bold',...
             'fontsize',13,'color',plot_c)
        txt_y=txt_y-0.11;

        ylabel('O_3 (DU)')

        %%% NO2 plot
        axes(fig_ax(2))
        
        ind=gbs_no2.year==yr;
        xx=mjd2k_to_date(gbs_no2.mjd2k);
        xx.Year=0;
        
        plot(xx(ind),gbs_no2.mean_vcd(ind),'s','color',plot_c,...
             'linewidth',1.2), hold on
            
        ylabel('NO_2 (molec/cm^2)')
        xlabel('Date (UTC)')

    end

    axes(fig_ax(1))
    grid on
    set(gca,'XTick',[61,70,80,92,101,111])
    title('GBS vertical columns from the PEARL Ridge Lab (Eureka, Canada, 80°N)')
    ylim([180,600])

    text(0.02,0.93,'Preliminary Results, Kristof Bognar, Ramina Alwarda, Kimberly Strong, UofT',...
         'units','normalized','fontsize',8,'fontweight','bold')

    axes(fig_ax(2))
    grid on
    set(gca,'XTick',[61,70,80,92,101,111])

    text(0.02,0.93,'Preliminary Results, Kristof Bognar, Ramina Alwarda, Kimberly Strong, UofT',...
         'units','normalized','fontsize',8,'fontweight','bold')

else
    
    figure
    set(gcf, 'Position', [100, 100, 1000, 500]);
    fig_ax = tight_subplot(1,1,[0.07,0.07],[0.125,0.1],[0.08,0.115]);

    txt_y=0.975;

    for yr=unique(gbs_o3.year)'


        ind=find(highlight==yr);

        if ~isempty(ind)
            plot_c=plot_colors(ind,:);
        else
            plot_c=[.5 .5 .5];
        end

        %%% oxone plot
        axes(fig_ax(1))

        ind=gbs_o3.year==yr;
        xx=mjd2k_to_date(gbs_o3.mjd2k);
        xx.Year=0;

        plot(xx(ind),gbs_o3.mean_vcd(ind),'s','color',plot_c,...
             'linewidth',1.2), hold on

        text(1.05,txt_y,num2str(yr),'units','normalized','fontweight','bold',...
             'fontsize',13,'color',plot_c)
        txt_y=txt_y-0.05;

        ylabel('Ozone column (Dobson units)')

        xlabel('Date (UTC)')

    end

    axes(fig_ax(1))
    grid on
    set(gca,'XTick',[61,70,80,92,101,111])
    title('Ozone total columns measured at the PEARL Ridge Lab (Eureka, Canada, 80°N)',...
          'fontsize',14)
    set(get(gca,'title'),'Position',[84 608 0])
    ylim([170,600])


    text(0.02,0.94,'Preliminary Results, Kristof Bognar, Ramina Alwarda, Kimberly Strong, UofT',...
         'units','normalized','fontsize',8,'fontweight','bold')

    set(gca,'XTick',[61,70,80,92,101,111])

    datetick('x', 'mmmm dd','keeplimits');

end
