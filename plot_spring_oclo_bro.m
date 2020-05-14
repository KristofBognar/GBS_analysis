%%% temporary code to get some OClO and BrO plots for a PAHA/ECCC telecon
%%% the datasets used are old and/or incomplete


highlight=[2000,2005,2007,2011,2014,2015,2020];

plot_colors=flipud(parula(length(highlight)));
plot_colors(highlight==2011,:)=plot_colors(highlight==2020,:);
plot_colors(highlight==2020,:)=[1 .1 0];

plot_gray=[.65 .65 .65];  

x_lim=[50,112];

if 1
    %% load data
    % quickly reprocessed 2011 data, with updated QDOAS project (includes BrO)
    load('/home/kristof/work/GBS/PEARL-GBS/2011_tmp_for_telecon/avg_twilight_for_oclo/oclo_2011.mat')
    data_new=data;
    data_new(1:220,:)=[];
    data_new(data_new.OClORMS>0.01,:)=[];

    % new 2020 data with updated QDOAS project (includes BrO)
    load('/home/kristof/work/GBS/PEARL-GBS/2020/avg_twilight_for_oclo/oclo_2020.mat')
    data(data.OClORMS>0.01,:)=[];
    data_new=[data_new;data];

    % old PGBS data from Cristen/me
    data_old_v1=[];
    load(['/home/kristof/work/GBS/PEARL-GBS/pre-2011_oclo/OClO_Eureka_2007_UToronto.mat']);
    data_old_v1=[data_old_v1;data];
    load(['/home/kristof/work/GBS/PEARL-GBS/pre-2011_oclo/OClO_Eureka_2010_UToronto.mat']);
    data=data(data.Year==2010,:);
    data_old_v1=[data_old_v1;data];

    for yr=2012:2016
        load(['/home/kristof/atmosp_servers/net/aurora/ground/eureka/gbs/pearl-gbs/'...
             num2str(yr) '/avg_twilight_for_oclo/OClO_Eureka_' num2str(yr) '_UToronto.mat']);

        data_old_v1=[data_old_v1;data];
    end

    data_old_v2=[];
    for yr=2017:2019
        load(['/home/kristof/atmosp_servers/net/aurora/ground/eureka/gbs/pearl-gbs/'...
             num2str(yr) '/avg_twilight_for_oclo/OClO_Eureka_' num2str(yr) '_UToronto.mat']);

        data_old_v2=[data_old_v2;data];
    end
end


%% approx detection limits

ind=abs(data_new.SZA-90)<1.1;

detlim_oclo=sqrt((mean(data_new.OClOSlErroclo(ind)))^2+(std(data_new.OClOSlColoclo(ind)))^2);

detlim_bro=sqrt((mean(data_new.BrOSlErrbro(ind)))^2+(std(data_new.BrOSlColbro(ind)))^2);


%% plot results

figure
set(gcf, 'Position', [100, 100, 1000, 600]);
fig_ax = tight_subplot(2,1,[0.07,0.07],[0.1,0.08],[0.08,0.13]);

txt_y=0.97;

axes(fig_ax(1))
plot([0,1],[1,1],'k--'), hold on % for legend

for yr=[2007,2010:2020]

    ind=find(highlight==yr);

    if ~isempty(ind)
        plot_c=plot_colors(ind,:);
    else
        plot_c=plot_gray;
    end
    
    if yr<2017 && yr~=2011
        
        ind=(abs(data_old_v1.SZA-90)<0.3 & data_old_v1.Year==yr);

        xx=ft_to_date(data_old_v1.DOY-1,data_old_v1.Year);
        xx.Year=0;
        yy=data_old_v1.OClO_DSCD;
        
    elseif yr>=2017 && yr<2020
        
        ind=(abs(data_old_v2.SZA-90)<0.1 & data_old_v2.Year==yr);

        xx=ft_to_date(data_old_v2.Fractionalday-1,data_old_v2.Year);
        xx.Year=0;
        yy=data_old_v2.NO2SlColoclo;
        
    elseif yr==2011 || yr==2020
        
        ind=(abs(data_new.SZA-90)<0.1 & data_new.Year==yr);

        xx=ft_to_date(data_new.Fractionalday-1,data_new.Year);
        xx.Year=0;
        yy=data_new.OClOSlColoclo;
        
    end
    
    plot(xx(ind),yy(ind),'s','color',plot_c,'linewidth',1.2), hold on
    
    text(1.05,txt_y,num2str(yr),'units','normalized','fontweight','bold',...
         'fontsize',13,'color',plot_c)
    txt_y=txt_y-0.11;


%     ind=abs(data_2020.SZA-90)<0.1;
% 
%     xx=ft_to_date(data_2020.Fractionalday-1,data_2020.Year);
%     xx.Year=0;
% 
%     plot(xx(ind),data_2020.OClOSlColoclo(ind),'s','color','r',...
%          'linewidth',1.2), hold on
% 
%     text(1.05,txt_y,'2020','units','normalized','fontweight','bold',...
%          'fontsize',13,'color','r')
%     txt_y=txt_y-0.11;

end

ylim([-1,3]*1e14)
xlim(x_lim)
ylabel('OClO dSCD (molec/cm^2)')

plot([x_lim(1)-1,x_lim(2)+1],[detlim_oclo,detlim_oclo],'k--')

xlim_arr=[61,70,80,92,101,111,122]-1;
grid on
set(gca,'XTick',xlim_arr)
set(gca,'XTickLabel',cellstr(ft_to_date(xlim_arr,0),'MMM dd'))

title('GBS-UV data: OClO and BrO slant columns at SZA=90°')
text(0.02,0.93,'To be reprocessed',...
     'units','normalized','fontsize',8,'fontweight','bold')

legend('Detection limit','location','northeast')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


axes(fig_ax(2))

for yr=[2011,2020]
    
    plot_c=plot_colors(highlight==yr,:);
    
    ind=(abs(data_new.SZA-90)<0.1 & data_new.Year==yr);

    xx=ft_to_date(data_new.Fractionalday-1,data_new.Year);
    xx.Year=0;
    yy=data_new.OClOSlColoclo;
            
    plot(xx(ind),yy(ind),'s','color',plot_c,'linewidth',1.2), hold on
    
end

ylim([-1,3]*1e14)
xlim(x_lim)
ylabel('BrO dSCD (molec/cm^2)')

plot([x_lim(1)-1,x_lim(2)+1],[detlim_bro,detlim_bro],'k--')
grid on
set(gca,'XTick',xlim_arr)
set(gca,'XTickLabel',cellstr(ft_to_date(xlim_arr,0),'MMM dd'))

text(0.02,0.93,'To be reprocessed',...
     'units','normalized','fontsize',8,'fontweight','bold')

xlabel('Date (UTC)')

