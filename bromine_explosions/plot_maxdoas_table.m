% to plot maxdoas results for PEARL-GBS


% load('/home/kristof/aurora/ground/eureka/gbs/pearl-gbs/2016/output/maxdoas_d66-131.mat')
% load('/home/kristof/aurora/ground/eureka/gbs/pearl-gbs/2010/QDOAS_output/maxdoas_d84-151.mat')
load('/home/kristof/work/GBS/PEARL-GBS/2018/QDOAS_output/maxdoas_d64-151.mat')
% load('/home/kristof/work/GBS/PEARL-GBS/2019/QDOAS_output/maxdoas_d64-151.mat')
% load('/home/kristof/work/GBS/PEARL-GBS/2017/QDOAS_output/maxdoas_d64-94.mat')
% load('/home/kristof/aurora/ground/eureka/gbs/pearl-gbs/2017/QDOAS_output/maxdoas_d66-94.mat')

%% control variables
x_date=1;
plot_o4=1;

msize=16;

% subtract=60; % 60 for leap years since data is in day of year
% subtract=59; % 59 for non-leap years since data is in day of year
subtract=1; % convert to fractional date

day1 = 64; 
day2 = 152;


%% set up the plot

% elevation angles in the file (might change from year to year)
elevs=unique(maxdoas.Elevviewingangle);

% colors for plotting
colors={'r.','g.','b.','y.','m.','c.','k.','ko'};

% empty legend cell
legends=cell(1,length(elevs));

% initialize figure
figure(99)
% figure('Position', [100, 100, 1410, 1050])
set(gcf, 'Position', [100, 100, 900, 600]);

if plot_o4, ax_bro=subplot(2,1,1); end
box on
hold on;
grid on
% grid minor
% ax.GridAlpha=0.4;
% grid minor

% select plotting axis
if x_date % datetime
    x_arr=maxdoas.DateTime;
else % DOY (-1 for fractional date)
    x_arr=maxdoas.Fractionalday-subtract;
end

%% loop over elevations
for i=1:length(elevs)
   
    % find indices with current elevation
    ind=find(maxdoas.Fractionalday >= day1 & maxdoas.Fractionalday <= day2 &...
             maxdoas.Elevviewingangle==elevs(i));
         
    % plot
    if i==8, msize_tmp=4; else msize_tmp=msize; end
    plot(x_arr(ind),maxdoas.BrOSlColbro(ind),colors{i},'markersize', msize_tmp);     
        
    % create legend entry
    legends{i}=[num2str(elevs(i)) '\circ'];
    
end 

ylim([-0.5e14 8e14]);
% xlabel('Fractional day')
ylabel('BrO DSCD (mol/cm^2)')
legend(legends,'location','northeast','Orientation','horizontal')

%%%%%%%%%%%%%
if plot_o4
    
    ax_o4=subplot(2,1,2);
    box on
    hold on;
    grid on
    % grid minor
    % ax.GridAlpha=0.4;
    % grid minor

    for i=1:length(elevs)

        % find indices with current elevation
        ind=find(maxdoas.Fractionalday >= day1 & maxdoas.Fractionalday <= day2 &...
                 maxdoas.Elevviewingangle==elevs(i));

        % plot
        if i==8, msize_tmp=4; else msize_tmp=msize; end
        plot(x_arr(ind),maxdoas.O4SlColo4(ind),colors{i},'markersize', msize_tmp);     

        % create legend entry
        legends{i}=[num2str(elevs(i)) '\circ'];

    end 
    
    ylim([-1e3 10e3]);
    % xlabel('Fractional day')
    ylabel('O_4 DSCD (mol/cm^2)')
    legend(legends,'location','northeast','Orientation','horizontal')
    
end

try linkaxes([ax_bro,ax_o3],'x'); end
try linkaxes([ax_bro,ax_o4],'x'); end


% set font on plots
set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman') 
% 
% f=gcf; 
% figpos=getpixelposition(f); 
% resolution=get(0,'ScreenPixelsPerInch'); 
% set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
% path='/home/kristof/work/summer_school/poster/'; 
% name='BrO_march_19-22'; 
% print(f,fullfile(path,name),'-dpng','-r300','-opengl') %save file

