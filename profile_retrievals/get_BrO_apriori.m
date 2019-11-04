function apriori=get_BrO_apriori(times)
%apriori=get_BrO_apriori(times): takes time of BrO profile, and returns
%corresponding a priori used in the BrO retrieval (surface pptv and scale height)

    surf_ppt=NaN(size(times));
    scale_h=NaN(size(times));
    
    x=datestr(times,'mmdd');
    
    yr=times(1).Year-1;
    
    for i=1:length(times)
        
        if times(i).Year == yr+1,
            load(['/home/kristof/work/profile_retrievals/profile_results/profile_details/'...
                  'prof_info_' num2str(times(i).Year) '.mat']);
            yr=times(i).Year;
        end
        
        surf_ppt(i)=str2double(apriori_BrO{find_in_cell(apriori_BrO(:,3),x(i,:)),1});
        scale_h(i)=str2double(apriori_BrO{find_in_cell(apriori_BrO(:,3),x(i,:)),2});
    end
    
    apriori=table();
    apriori.surf_ppt=surf_ppt;
    apriori.scale_h=scale_h;
    
end


