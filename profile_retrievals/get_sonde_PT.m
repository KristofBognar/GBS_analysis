function out = get_sonde_PT( year, day_range, dmean, plot_prof )
%get_sonde_PT( year, day_range, dmean, plot_prof ) info about radiosonde PT profiles
%   Returns either daily mean temperature info for retrievals (dmean=1)
%       or profile info corresponding to individual sondes (dmean=0)
%   Individual profiles can be plotted (plot_prof=1) -- click on plot to
%       advance to next profile

% year=2018;
% day_range=[64,151];
% 
% plot_prof=0;

if nargin==1
   
    day_range=[0,366];
    dmean=0;
    plot_prof=0;
    
end

load(['/home/kristof/work/radiosonde/Eureka/radiosonde_' num2str(year) '_interp.mat']);

% launch times indicated as noon and midnight in GRAW files; 
% launches are probably 11 someting and 23 something -- subtract 30 min
% from fractional times
if year>=2019; x=x-0.01; end

times=ft_to_date(x,year); 

row_ind=1;

times_all=datetime;
dT_all=[];

if dmean
    %% Find RL to sea level temperature difference (mean value for each day)
    % Strong inversion (scale height = 0.5 km) when dT>7

    % one T difference for entire day
    for i=day_range(1):day_range(2)

        % find all sondes on given day
        inds=find(x>i-1 & x<=i); 

        if ~isempty(inds)
            % save day
            times_all(row_ind,1)=times(inds(1));
            % save mean T difference between ridge lab and sea level
            dT_all(row_ind,1)=mean(T_arr(60,inds)-T_arr(1,inds));
            row_ind=row_ind+1;
        end
    end
    
    out=table;
    out.date=times_all;
    out.date.Format='dd/MM/yyyy';
    out.dT=dT_all;

    out.scale_h=ones(size(out.dT));
    out.scale_h(out.dT>7)=0.5;
    
    
else
    %% Return info from individual profiles
   
    % T difference from individual profiles
    T_600=T_arr(60,:)';
    T_400=T_arr(40,:)';
    T_200=T_arr(20,:)';
    T_0=T_arr(1,:)';
    dT_all=T_600-T_0;

    P_600=P_arr(60,:)';
    P_0=P_arr(1,:)';
    
    times_all=times';
    
    dT_all(x<day_range(1) | x>day_range(2))=[];
    T_0(x<day_range(1) | x>day_range(2))=[];
    T_200(x<day_range(1) | x>day_range(2))=[];
    T_400(x<day_range(1) | x>day_range(2))=[];
    T_600(x<day_range(1) | x>day_range(2))=[];
    
    P_0(x<day_range(1) | x>day_range(2))=[];
    P_600(x<day_range(1) | x>day_range(2))=[];
    
    times_all(x<day_range(1) | x>day_range(2))=[];
    
    out=table;
    out.date=times_all;
    out.date.Format='dd/MM/yyyy HH:mm';
    out.T_0=T_0-273.15;
    out.T_200=T_200-273.15;
    out.T_400=T_400-273.15;
    out.T_600=T_600-273.15;
    out.dT=dT_all;

    out.P_0=P_0;
    out.P_600=P_600;
end

%% Plot all temperature profiles

if plot_prof
    
    inds=find(x>day_range(1) & x<day_range(2)+1);

    times=ft_to_date(x(inds),year);

    x_lim=[-40,-10];
    
    for i=1:length(inds)

        % ridge lab altitude
        plot([-60,60],[0.61,0.61],'k--'), hold on

        % redraw previous line in gray
        if i>1
            delete(h)
            plot(T_arr(:,inds(i-1))-273.15,y/1000,'color',[0,0,0]+0.8,'linewidth',2), hold on
        end

        h=plot(T_arr(:,inds(i))-273.15,y/1000,'b-','linewidth',2);
        ylim([0,2])
        
        % make x limits adaptive
        if max(T_arr(1:200,inds(i)))-273.15>x_lim(2)
            x_lim(2)=max(T_arr(1:200,inds(i)))-272.15;
        end
        if min(T_arr(1:200,inds(i)))-273.15<x_lim(1)
            x_lim(1)=min(T_arr(1:200,inds(i)))-274.15;
        end
        xlim(x_lim)


        title([datestr(times(i),'mmm dd - HH:SS') ', day ' num2str(floor(x(inds(i)))+1)])

        xlabel('Temperature (°C)')
        ylabel('Altitude (km)')
        
        % pause loop until figure is clicked on
        try
            tmp=1;
            while tmp % loop so key presses (return 1) are not accepted
                tmp=waitforbuttonpress;
            end
        catch
            % returns error if figure is closed, exit when that happens
            return
        end
    end

end

% end





