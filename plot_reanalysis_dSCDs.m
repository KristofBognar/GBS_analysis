
uv=false;
vis=false;
pvis=true;

if vis 
    load('o3_X_U_1999-2016.mat') % all saved variables are sorted by measurement time
    load('full_U_1999-2016.mat')
    load('full_U_2008-2010_UV.mat')

    %% plot total ozone RMS, with and without X.xs

    % filter large RMS and SZA + select year
    ind=find(data(:,146)<0.01 & data(:,7)<91 & data(:,2)==2011);
    ind_o3x=find(o3x_rms<0.01 & o3x_data(:,7)<91 & o3x_data(:,2)==2011);

    figure(1)
    plot(data(ind,3),data(ind,146)), hold on
    plot(o3x_data(ind_o3x,3),o3x_rms(ind_o3x),'r--'), hold on
    plot(UTGBS2011fixedref0622(:,3),UTGBS2011fixedref0622(:,196),'k'), hold on
    % plot(mjd(ind),data(ind,146)), hold on
    % plot(o3x_mjd(ind_o3x),o3x_rms(ind_o3x),'r--'), hold on

    figure(2)
    % plot(mjd(ind),data(ind,157),'bo','markersize',4), hold on
    % plot(o3x_mjd(ind_o3x),o3x_dscd(ind_o3x),'r.','markersize',10), hold on
    plot(data(ind,3),data(ind,157),'bo','markersize',4), hold on
    plot(o3x_data(ind_o3x,3),o3x_dscd(ind_o3x),'r.','markersize',10), hold on
    plot(UTGBS2011fixedref0622(:,3),UTGBS2011fixedref0622(:,207),'k.','markersize',10), hold on

    %% plot yearly ozone DSCDs and RMS, with and without X.xs

    % year=[];
    % mean_dscd=[];
    % mean_rms=[];
    % mag_X=[];
    % 
    % ind_spring=find(data(:,3)>=0 & data(:,3)<=120);
    % ind_spring_o3x=find(o3x_data(:,3)>=0 & o3x_data(:,3)<=120);
    % 
    % for i=1999:2016
    %     
    %    ind_year=find(data(:,2)==i & data(:,146)<0.01 & data(:,7)<91);
    %    ind_year_o3x=find(o3x_data(:,2)==i & o3x_rms<0.01 & o3x_data(:,7)<91);
    % 
    % %    ind_year=intersect(ind_year,ind_spring);
    % %    ind_year_o3x=intersect(ind_year_o3x,ind_spring_o3x);
    %    
    %    if isempty(ind_year), continue, end
    %    
    %    year=[year; i];
    %    mean_dscd=[mean_dscd; [nanmean(data(ind_year,157)), nanmean(o3x_dscd(ind_year_o3x))] ];
    %    mean_rms=[mean_rms; [nanmean(data(ind_year,146)), nanmean(o3x_rms(ind_year_o3x))] ];
    %    mag_X=[mag_X; nanmean(abs(o3x_data(ind_year_o3x,27)))];
    %    
    % end
    % 
    % figure(2)
    % 
    % subplot(2,3,[3,6])
    % plot(year,mean_dscd(:,1),'bo','linewidth',1.2,'markersize',8), hold on
    % plot(year,mean_dscd(:,2),'rx','linewidth',1.2,'markersize',8)
    % xlim([1998,2017])
    % ylim([2,8]*1e19)
    % xlabel('Year')
    % ylabel('Mean O_3 DSCD')
    % legend('No X.xs','With X.xs','location','northwest')
    % 
    % subplot(2,3,[1,2])
    % plot(year,mean_rms(:,1),'bo','linewidth',1.2,'markersize',8), hold on
    % plot(year,mean_rms(:,2),'rx','linewidth',1.2,'markersize',8)
    % xlim([1998,2017])
    % ylabel('Mean O_3 RMS')
    % legend('No X.xs','With X.xs','location','northwest')
    % 
    % subplot(2,3,[4,5])
    % plot(year,mag_X,'kx--','linewidth',1.2)
    % xlim([1998,2017])
    % xlabel('Year')
    % ylabel('Mean X DSCD')

    %% plot yearly NO2 RMS

    % year=[];
    % mean_dscd=[];
    % mean_rms=[];
    % 
    % ind_spring=find(data(:,3)>=0 & data(:,3)<=120);
    % ind_spring_uv=find(uv_data(:,3)>=0 & uv_data(:,3)<=120);
    % 
    % for i=1999:2016
    %     
    %    ind_year=find(data(:,2)==i & data(:,170)<0.01 & data(:,7)<91);
    %    ind_year_uv=find(uv_data(:,2)==i & uv_data(:,14)<0.01 & uv_data(:,7)<91);
    % 
    %    ind_year=intersect(ind_year,ind_spring);
    %    ind_year_uv=intersect(ind_year_uv,ind_spring_uv);
    %    
    %    if isempty(ind_year), continue, end
    %    
    %    year=[year; i];
    %    mean_dscd=[mean_dscd; [nanmean(data(ind_year,179)),nanmean(uv_data(ind_year_uv,21))]];
    %    mean_rms=[mean_rms; [nanmean(data(ind_year,170)),nanmean(uv_data(ind_year_uv,14))]];
    %    
    % end
    % 
    % figure(2)
    % 
    % subplot(222)
    % plot(year,mean_dscd(:,1),'bo','linewidth',1.2,'markersize',8), hold on
    % plot(year,mean_dscd(:,2),'rx','linewidth',1.2,'markersize',8), hold on
    % xlim([1998,2017])
    % ylabel('Mean NO_2 DSCD')
    % legend('VIS','UV')
    % % title('Spring average')
    % 
    % subplot(224)
    % plot(year,mean_rms(:,1),'bo','linewidth',1.2,'markersize',8), hold on
    % plot(year,mean_rms(:,2),'rx','linewidth',1.2,'markersize',8), hold on
    % xlim([1998,2017])
    % xlabel('Year')
    % ylabel('Mean NO_2 RMS')
    % legend('VIS','UV')
    % 
    
elseif uv || pvis
    
    if uv
        load part_P_UV.mat
        if false 
            % ozone
            rms=14;
            dscd=21;
        else
            % no2
            rms=42;
            dscd=49;
        end
    
    elseif pvis
        load part_P_VIS.mat
        if true 
%             % ozone
%             rms=146;
%             dscd=157;
            % ozone with X.xs
            rms=170;
            dscd=181;
        else
            % no2
            rms=196;
            dscd=205;
        end
    end
    
    sza_range=[8,91];
    
    
    
    ind07=find(pgbs2007(:,rms)<0.05 & pgbs2007(:,7)<sza_range(2) & ...
               pgbs2007(:,7)>sza_range(1));
    ind08=find(pgbs2008(:,rms)<0.05 & pgbs2008(:,7)<sza_range(2) & ...
               pgbs2008(:,7)>sza_range(1));
    ind09=find(pgbs2009(:,rms)<0.05 & pgbs2009(:,7)<sza_range(2) & ...
               pgbs2009(:,7)>sza_range(1));
    ind10=find(pgbs2010(:,rms)<0.05 & pgbs2010(:,7)<sza_range(2) & ...
               pgbs2010(:,7)>sza_range(1));
%     ind10_3=find(pgbs2010_3spec(:,rms)<0.05 & pgbs2010_3spec(:,7)<sza_range(2)...
%                  & pgbs2010_3spec(:,7)>sza_range(1));
%     ind10_5=find(pgbs2010_5spec(:,rms)<0.05 & pgbs2010_5spec(:,7)<sza_range(2)...
%                  & pgbs2010_5spec(:,7)>sza_range(1));
%     ind11=find(pgbs2011(:,rms)<0.05 & pgbs2011(:,7)<sza_range(2) & ...
%                pgbs2011(:,7)>sza_range(1));
%     ind12=find(pgbs2012(:,rms)<0.05 & pgbs2012(:,7)<sza_range(2) & ...
%                pgbs2012(:,7)>sza_range(1));
%     ind13=find(pgbs2013(:,rms)<0.05 & pgbs2013(:,7)<sza_range(2) & ...
%                pgbs2013(:,7)>sza_range(1));
%     ind14=find(pgbs2014(:,rms)<0.05 & pgbs2014(:,7)<sza_range(2) & ...
%                pgbs2014(:,7)>sza_range(1));
%     ind15=find(pgbs2015(:,rms)<0.05 & pgbs2015(:,7)<sza_range(2) & ...
%                pgbs2015(:,7)>sza_range(1));
    
    plot(pgbs2007(ind07,3),pgbs2007(ind07,rms),'kx'), hold on

end


    
    