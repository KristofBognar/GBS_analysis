function prof_len=match_prof_length(times)
%prof_len=match_prof_length(times): takes time of BrO profile, and returns
%corresponding time interval (total time of the corresponding elevation
%scan)

    prof_len=NaN(size(times));
    x=datestr(times,'mmdd');
    
    yr=times(1).Year-1;
    
    for i=1:length(times)
        
        if times(i).Year == yr+1,
            load(['/home/kristof/work/profile_retrievals/profile_results/profile_details/'...
                  'prof_info_' num2str(times(i).Year) '.mat']);
            yr=times(i).Year;
        end
        
        prof_len(i)=str2double(daily_times{find_in_cell(daily_times(:,4),x(i,:)),3});
    end
    
end


