% compare VCDs across the entire period to Cristen's reanalysis

tg=1;

load('/home/kristof/work/GBS/VCD_results/VCDs_Cristen_JUNE2011.mat');

if tg==1

    load('/home/kristof/work/GBS/VCD_results/UT-GBS_O3_VCD_all.mat');
    utgbs=reanalysis;

    load('/home/kristof/work/GBS/VCD_results/PEARL-GBS_O3_VCD_all.mat');
    pgbs=reanalysis;
    
    %% find coincidences
    
    ut_diff=NaN(length(utgbs.year));
    p_diff=NaN(length(pgbs.year));
    
    time_u=utgbs.year+((utgbs.fd-1)./365);
%     time_uc=vcd_u1_o3_ndacc_filt(:,1)+((vcd_u1_o3_ndacc_filt(:,4)-1)./365);

    time_p=pgbs.year+((pgbs.fd-1)./365);
%     time_pc=vcd_p1_o3_ndacc_filt(:,1)+((vcd_p1_o3_ndacc_filt(:,4)-1)./365);
    
    for i=1:length(utgbs.year)

        ind=find(vcd_u1_o3_ndacc_filt(:,1)==utgbs.year(i) & ...
                 vcd_u1_o3_ndacc_filt(:,2)==utgbs.day(i) & ...
                 vcd_u1_o3_ndacc_filt(:,3)==utgbs.ampm(i));
         
        if ~isempty(ind)     
            ut_diff(i)=(utgbs.mean_vcd(i)/ vcd_u1_o3_ndacc_filt(ind,13)) -1;
        end
        
    end

    for i=1:length(pgbs.year)

        ind=find(vcd_p1_o3_ndacc_filt(:,1)==pgbs.year(i) & ...
                 vcd_p1_o3_ndacc_filt(:,2)==pgbs.day(i) & ...
                 vcd_p1_o3_ndacc_filt(:,3)==pgbs.ampm(i));
             
        if ~isempty(ind)                  
            p_diff(i)=(pgbs.mean_vcd(i)/ vcd_p1_o3_ndacc_filt(ind,13)) -1;
        end

    end
    
    figure(1)
        
    subplot(211)    
    plot(time_u,ut_diff*100, 'k.')
    title('Ozone relative differences compared to Cristens reanalysis')
    ylabel('UT-GBS rel diff (%)')
        
    subplot(212)
    plot(time_p,p_diff*100, 'k.')
    ylabel('PEARL-GBS rel diff (%)')    
    xlabel('Fractional year (UTC)')
    
    
    
elseif tg==2

    load('/home/kristof/work/GBS/VCD_results/UT-GBS_NO2_VCD_all.mat');
    utgbs=reanalysis;

    load('/home/kristof/work/GBS/VCD_results/PEARL-GBS_NO2_VCD_all.mat');
    pgbs=reanalysis;
    
    %% find coincidences
    
    ut_diff=NaN(length(utgbs.year));
    p_diff=NaN(length(pgbs.year));
    
    time_u=utgbs.year+((utgbs.fd-1)./365);
%     time_uc=vcd_u1_o3_ndacc_filt(:,1)+((vcd_u1_o3_ndacc_filt(:,4)-1)./365);

    time_p=pgbs.year+((pgbs.fd-1)./365);
%     time_pc=vcd_p1_o3_ndacc_filt(:,1)+((vcd_p1_o3_ndacc_filt(:,4)-1)./365);
    
    for i=1:length(utgbs.year)

        ind=find(vcd_u1_no2_ndacc_filt(:,1)==utgbs.year(i) & ...
                 vcd_u1_no2_ndacc_filt(:,2)==utgbs.day(i) & ...
                 vcd_u1_no2_ndacc_filt(:,3)==utgbs.ampm(i));
         
        if ~isempty(ind)     
            ut_diff(i)=(utgbs.mean_vcd(i)/ vcd_u1_no2_ndacc_filt(ind,13)) -1;
        end
        
    end

    for i=1:length(pgbs.year)

        ind=find(vcd_p1_no2_ndacc_filt(:,1)==pgbs.year(i) & ...
                 vcd_p1_no2_ndacc_filt(:,2)==pgbs.day(i) & ...
                 vcd_p1_no2_ndacc_filt(:,3)==pgbs.ampm(i));
             
        if ~isempty(ind)                  
            p_diff(i)=(pgbs.mean_vcd(i)/ vcd_p1_no2_ndacc_filt(ind,13)) -1;
        end

    end
    
    figure(1)
    
    subplot(211)
    plot(time_u,ut_diff*100, 'k.')
    title('NO_2 relative differences compared to Cristens reanalysis')
    ylabel('UT-GBS rel diff (%)')
    
    subplot(212)
    plot(time_p,p_diff*100, 'k.')
    ylabel('PEARL-GBS rel diff (%)')    
    xlabel('Fractional year (UTC)')
    
end


