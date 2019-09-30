% plot measured and modelled dscds for gived day


allelevs=false;

figure(1)
subplot(211)

plot(-1,1,'k','marker','square','linestyle','none'), hold on
plot(-1,1,'k','marker','x','linestyle','none'), hold on

if allelevs
    
    load 20160320_noadaptive-1.mat
    
    elevs=[-1,0,1,2,5,10,15,30];
    color=['r','g','b','y','m','c','k','k'];

    for i=1:length(elevs)
        plot(-1, elevs(i) ,'marker','none','linestyle','none'), hold on
    end
    legend('Meas','Model','\color{red}-1','\color{green}0','\color{blue}1',...
        '\color{yellow}2','\color{magenta}5',...
        '\color{cyan}10','\color{black}15','\color{black}30',...
        'location','eastoutside');
    
else
    
    load 20160320_noadaptive_no0_1iter.mat

    elevs=[1,2,5,10,15,30];
    color=['b','y','m','c','k','k'];
    
    for i=1:length(elevs)
        plot(-1, elevs(i) ,'marker','none','linestyle','none'), hold on
    end
    legend('Meas','Model','\color{blue}1','\color{yellow}2','\color{magenta}5',...
        '\color{cyan}10','\color{black}15','\color{black}30',...
        'location','eastoutside');
end

for i=1:length(elevs)
    
    ind=find(dscd(:,2)==elevs(i));
    
    time=24*(ft_dscd(ind)-79);
    
    plot(time,dscd(ind,5),color(i),...
        'marker','square','linestyle','none'), hold on
    plot(time,dscd(ind,7),color(i),...
        'marker','x','linestyle','none'), hold on

    xlim([13,23])
end

xlabel('Hours of March 20th, 2016 (UTC)')
ylabel('BrO DSCD (molec/cm^2)')
