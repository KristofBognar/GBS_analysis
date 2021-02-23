

instr='UT-GBS';
% instr='PEARL-GBS';

tg='O3';
% tg='NO2';
% tg='NO2_UV';

switch instr
    case 'UT-GBS'
        if ismac
            load(['/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/VCD_results/UT-GBS_' tg '_VCD_all.mat'])
        elseif isunix
            load(['/home/kristof/work/GBS/VCD_results/UT-GBS_' tg '_VCD_all.mat'])
        end
    case 'PEARL-GBS'
        if ismac
            load(['/Users/raminaalwarda/Downloads/PEARL-GBS_' tg '_VCD_all.mat'])
        elseif isunix
%         load(['/home/kristof/work/GBS/VCD_results/VCD_with_bad_rcd/PEARL-GBS_' tg '_VCD_all.mat'])
            load(['/home/kristof/work/GBS/VCD_results/PEARL-GBS_' tg '_VCD_all.mat'])
        end    
end


years=unique(reanalysis.year);

colors=jet(length(years));
% colors=flipud(hsv(length(years)));

legendstr={};

figure()
% subplot(223)

for i=1:length(years)

    ind=find(reanalysis.year==years(i));
    scatter(reanalysis.fd(ind)-1,reanalysis.mean_vcd(ind),20,colors(i,:),...
            'filled', 'MarkerFaceColor',colors(i,:)), hold on

    legendstr{i}=num2str(years(i));

end

legend(legendstr,'location','eastoutside')

xlabel('Fractional day (UTC)')

if ~strcmp(tg,'NO2_UV')
    ylabel([tg(1:end-1) '_' tg(end)  ' VCD (molec/cm^2)'])
else
    ylabel('NO_2 VCD (molec/cm^2)')
end

f=gcf; 
f.Units = 'pixels';

set(f, 'Position', [100, 100, 900, 450]); 
set(findall(gcf,'-property','FontSize'),'FontSize',13.0) 
set(findall(gcf,'-property','FontName'),'FontName','Arial') 
