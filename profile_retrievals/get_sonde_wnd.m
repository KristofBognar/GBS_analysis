function out = get_sonde_wnd( year, day_range, plot_prof )
%get_sonde_wnd( year, day_range, dmean, plot_prof ) Get radiosonde wind info

% year=2018;
% day_range=[64,151];
% 
% plot_prof=0;

if nargin==1
   
    day_range=[0,366];
    plot_prof=0;
    
end

load(['/home/kristof/work/radiosonde/Eureka/radiosonde_wnd_' num2str(year) '_interp.mat']);

times=ft_to_date(x,year);

%% Return info from individual profiles

% T difference from individual profiles
wspd_600=wspd_arr(60,:)';
wspd_400=wspd_arr(40,:)';
wspd_200=wspd_arr(20,:)';
wspd_0=wspd_arr(1,:)';

wdir_600=wdir_arr(60,:)';
wdir_400=wdir_arr(40,:)';
wdir_200=wdir_arr(20,:)';
wdir_0=wdir_arr(1,:)';

times_all=times';

wspd_0(x<day_range(1) | x>day_range(2))=[];
wspd_200(x<day_range(1) | x>day_range(2))=[];
wspd_400(x<day_range(1) | x>day_range(2))=[];
wspd_600(x<day_range(1) | x>day_range(2))=[];

wdir_0(x<day_range(1) | x>day_range(2))=[];
wdir_200(x<day_range(1) | x>day_range(2))=[];
wdir_400(x<day_range(1) | x>day_range(2))=[];
wdir_600(x<day_range(1) | x>day_range(2))=[];

times_all(x<day_range(1) | x>day_range(2))=[];

out=table;
out.date=times_all;
out.date.Format='dd/MM/yyyy HH:mm';

out.wspd_0=wspd_0;
out.wspd_200=wspd_200;
out.wspd_400=wspd_400;
out.wspd_600=wspd_600;

out.wdir_0=wdir_0;
out.wdir_200=wdir_200;
out.wdir_400=wdir_400;
out.wdir_600=wdir_600;


%% Plot all temperature profiles

if plot_prof
    
    inds=find(x>day_range(1) & x<day_range(2)+1);

    times=ft_to_date(x(inds),year);

    for i=1:length(inds)

        % ridge lab altitude
        plot([-1,22],[0.61,0.61],'k--'), hold on

        % redraw previous line in gray
        if i>1
            delete(h)
            plot(wspd_arr(:,inds(i-1)),y/1000,'color',[0,0,0]+0.8,'linewidth',2), hold on
        end

        h=plot(wspd_arr(:,inds(i)),y/1000,'b-','linewidth',2);
        ylim([0,2])
        xlim([-1,22])

        title([datestr(times(i),'mmm dd - HH:SS') ', day ' num2str(floor(x(inds(i)))+1)])

        xlabel('Wind speed (m/s)')
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





