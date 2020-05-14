%%% plot yearly ozone and NO2 in the spring, to look at ozone depletion
%%% VCD plot modelled on Adams et al., 2012

du=2.687e16;

curr_yr=2020;

vcd_plot=1; % 1: plot ozone and NO2, 2: ozone only
no2_diffs=0;

instr='UT-GBS';
% instr='PEARL-GBS';

highlight=[2000,2005,2007,2011,2014,2015,2020];
x_lim=[50,122];

% plot_colors=flipud([[1 0 0];...
%              [0 0 .6];...
%              [0 0 .3];...
%              [0 0 1];...
%              [0 1 0];...
%              [1 1 0];...
%              [1 0 1];...
%                     ]);

plot_gray=[.65 .65 .65];  

diurnal_bins=55:5:120;

plot_colors=flipud(parula(length(highlight)));
plot_colors(highlight==2011,:)=plot_colors(highlight==2020,:);
plot_colors(highlight==2020,:)=[1 .1 0];

% merged GBS data, up to 2017 (from merged VCDs, already filtered)
load('/home/kristof/work/satellite_validation/GBS_O3.mat')
load('/home/kristof/work/satellite_validation/GBS_NO2.mat')
load('/home/kristof/work/satellite_validation/GBS_NO2_UV.mat')

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


%% add later NO2-UV data
% load merged files only, since those are filtered
load('/home/kristof/work/GBS/VCD_results/PEARL-GBS_NO2_UV_VCD_all.mat')
gbs_no2uv=[gbs_no2uv; reanalysis(reanalysis.year>2017,:)];

% check if some data is still RD (no standard processing yet)
check=max(reanalysis.year)+1:curr_yr;

for yr=check
    
    load(['/home/kristof/work/GBS/VCD_results/NDACC_RD/PEARL-GBS_NO2_UV_VCD_' num2str(yr) 'all.mat'])
    gbs_no2uv=[gbs_no2uv; reanalysis];
    
end


%% clean up data

% cut off at specified date
gbs_o3(gbs_o3.day>x_lim(2),:)=[];
gbs_no2(gbs_no2.day>x_lim(2),:)=[];
gbs_no2uv(gbs_no2uv.day>x_lim(2),:)=[];

% remove single twilights for NO2 (should be very few, if any)
yd=gbs_no2.year*1000+gbs_no2.day;
tmp=unique(yd);

count=histc(yd,tmp);
tmp=tmp(count==1);

if ~isempty(tmp)
    for i=tmp'
        gbs_no2(yd==i,:)=[];
        yd(yd==i,:)=[];
    end
end

% same for NO2-UV
yd=gbs_no2uv.year*1000+gbs_no2uv.day;
tmp=unique(yd);

count=histc(yd,tmp);
tmp=tmp(count==1);

if ~isempty(tmp)
    for i=tmp'
        gbs_no2uv(yd==i,:)=[];
        yd(yd==i,:)=[];
    end
end

%% plot results

if vcd_plot==1

    figure
    set(gcf, 'Position', [100, 100, 1000, 600]);
    fig_ax = tight_subplot(2,1,[0.07,0.07],[0.1,0.08],[0.08,0.13]);

    txt_y=0.97;

    for yr=unique(gbs_o3.year)'


        ind=find(highlight==yr);

        if ~isempty(ind)
            plot_c=plot_colors(ind,:);
        else
            plot_c=plot_gray;
        end

        %%% ozone plot
        axes(fig_ax(1))

        ind=gbs_o3.year==yr;
        xx=mjd2k_to_date(gbs_o3.mjd2k);
        xx.Year=0;

        plot(xx(ind),gbs_o3.mean_vcd(ind),'s','color',plot_c,...
             'linewidth',1.2), hold on

        text(1.05,txt_y,num2str(yr),'units','normalized','fontweight','bold',...
             'fontsize',13,'color',plot_c)
        txt_y=txt_y-0.11;

        ylabel('O_3 VCD (DU)')

        %%% NO2 plot
        axes(fig_ax(2))
        
        ind=gbs_no2.year==yr;
        xx=mjd2k_to_date(gbs_no2.mjd2k);
        xx.Year=0;
        
        plot(xx(ind),gbs_no2.mean_vcd(ind),'s','color',plot_c,...
             'linewidth',1.2), hold on
            
        ylabel('NO_2 VCD (molec/cm^2)')
        xlabel('Date (UTC)')

    end

    axes(fig_ax(1))
    grid on
    set(gca,'XTick',[61,70,80,92,101,111,122])
%     title('GBS vertical columns from the PEARL Ridge Lab (Eureka, Canada, 80°N)')
    title('GBS-Vis data: vertical columns of O_3 and NO_2')
    ylim([180,600])

%     text(0.02,0.93,'Preliminary Results, Kristof Bognar, Ramina Alwarda, Kimberly Strong, UofT',...
%          'units','normalized','fontsize',8,'fontweight','bold')

    axes(fig_ax(2))
    grid on
    set(gca,'XTick',[61,70,80,92,101,111,122])

%     text(0.02,0.93,'Preliminary Results, Kristof Bognar, Ramina Alwarda, Kimberly Strong, UofT',...
%          'units','normalized','fontsize',8,'fontweight','bold')

elseif vcd_plot==2
    
    figure
    set(gcf, 'Position', [100, 100, 1000, 500]);
    fig_ax = tight_subplot(1,1,[0.07,0.07],[0.125,0.1],[0.08,0.115]);

    txt_y=0.975;

    for yr=unique(gbs_o3.year)'


        ind=find(highlight==yr);

        if ~isempty(ind)
            plot_c=plot_colors(ind,:);
        else
            plot_c=plot_gray;
        end

        %%% ozone plot
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

if no2_diffs
    
    
    figure
    set(gcf, 'Position', [100, 100, 1000, 600]);
    fig_ax = tight_subplot(2,1,[0.07,0.07],[0.1,0.08],[0.08,0.13]);

    txt_y=0.97;

    %%% plot diurnal diffs for both UV and vis
    gray_yr=setdiff(unique(gbs_no2.year)',highlight);

    mm_max=NaN(1,length(diurnal_bins));
    mm_min=mm_max;
    xx=diurnal_bins-0.25;
%     xx=ft_to_date(xx,yr);
%     xx.Year=0;

    for uvvis=1:2
        
        if uvvis==1
            no2_table=gbs_no2;
            axes(fig_ax(1))
        else
            no2_table=gbs_no2uv; % starts in 2007
            axes(fig_ax(2))
        end
        
%         for yr=gray_yr
% 
%             %%% diurnal differences plot, NO2-VIS
%             ind_am=(no2_table.year==yr & no2_table.ampm==0);
%             ind_pm=(no2_table.year==yr & no2_table.ampm==1);
% 
%             yy=no2_table.mean_vcd(ind_pm)-no2_table.mean_vcd(ind_am);
% 
%             % mean diurnal diffs
%             for i=1:length(diurnal_bins)
%                 ind=(abs(no2_table.day(ind_am)-diurnal_bins(i))<3); % +- 2 days = 5 day bins
%                 if ~isempty(ind) 
%                     mm_max(i)=max([mean(yy(ind)),mm_max(i)]);
%                     mm_min(i)=min([mean(yy(ind)),mm_min(i)]);
%                 end
%             end
% 
%         end
% 
%         fill([xx,fliplr(xx)],[mm_max,fliplr(mm_min)],plot_gray+0.2,'LineStyle','none'), hold on

            
%         for yr=highlight
        for yr=unique(gbs_no2.year)'

%             plot_c=plot_colors(highlight==yr,:);
            ind=find(highlight==yr);

            if ~isempty(ind)
                plot_c=plot_colors(ind,:);
                ls='s-';
            else
                plot_c=plot_gray;
                ls='s';
            end

            %%% diurnal differences plot, NO2-VIS
            ind_am=(no2_table.year==yr & no2_table.ampm==0);
            ind_pm=(no2_table.year==yr & no2_table.ampm==1);

            yy=no2_table.mean_vcd(ind_pm)-no2_table.mean_vcd(ind_am);

            mm=NaN(1,length(diurnal_bins));
            % mean diurnal diffs
            for i=1:length(diurnal_bins)
                ind=(abs(no2_table.day(ind_am)-diurnal_bins(i))<3); % +- 2 days = 5 day bins
                if ~isempty(ind) 
                    mm(i)=mean(yy(ind));
                end
            end

            plot(xx,mm,ls,'color',plot_c,'linewidth',1.2,'markerfacecolor',plot_c), hold on

        end
        
        ylim([-3,12]*1e14)
        xlim(x_lim)
        
    end
    
    xlim_arr=[61,70,80,92,101,111,122]-1;
    
    axes(fig_ax(1))
    grid on
    set(gca,'XTick',xlim_arr)
    set(gca,'XTickLabel',cellstr(ft_to_date(xlim_arr,0),'MMM dd'))
    title('GBS Vis and UV data: diurnal variation of NO_2 (sunset minus sunrise, 5-day mean)')
    ylabel('\DeltaNO_2 (molec/cm^2)')

%     text(0.02,0.93,'Visible measurements (NDACC)',...
%          'units','normalized','fontsize',8,'fontweight','bold')

    %%% plot years on the side
    axes(fig_ax(1))
    for yr=unique(gbs_no2.year)'

        ind=find(highlight==yr);

        if ~isempty(ind)
            plot_c=plot_colors(ind,:);
        else
            plot_c=plot_gray;
        end

        text(1.05,txt_y,num2str(yr),'units','normalized','fontweight','bold',...
             'fontsize',13,'color',plot_c)
        txt_y=txt_y-0.11;
         
    end
    
%     legend(['Full range of diurnal variability for' char(10) 'years with no vortex measurements'])

     
    axes(fig_ax(2))
    grid on
    set(gca,'XTick',xlim_arr)
    set(gca,'XTickLabel',cellstr(ft_to_date(xlim_arr,0),'MMM dd'))
    text(0.02,0.93,'UV measurements available for 2007-2020',...
         'units','normalized','fontsize',8,'fontweight','bold')
    
    xlabel('Date (UTC)')
    ylabel('\DeltaNO_2-UV (molec/cm^2)')
     
end




