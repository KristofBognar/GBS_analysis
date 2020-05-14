% find measurements within the polar vortex using DMPs along GBS LoS

if 1
    % use ozone data, will have the most datapoints
    % merged GBS data, up to 2017 (from merged VCDs, already filtered)
    load('/home/kristof/work/satellite_validation/GBS_O3.mat')
    
    gbs_o3(gbs_o3.day>122,:)=[];
    
    disp('Assumin that DMPS are available up to 2017 only')

    spv = match_DMP_DOAS( gbs_o3.fractional_time, gbs_o3.year, 'O3_VIS', 19.5 );
end

figure

legend_cell=[];
cnt=1;
for yr=unique(gbs_o3.year)'
    
    ind=gbs_o3.year==yr;
    if ~isempty(ind)
        plot(gbs_o3.fractional_time(ind),spv(ind),'s-','linewidth',1.2), hold on
        legend_cell{cnt}=num2str(yr);
        cnt=cnt+1;
    end
end

plot([50,123],[1.6,1.6]*1e-4,'k--')
plot([50,123],[1.2,1.2]*1e-4,'k--')
legend(legend_cell)