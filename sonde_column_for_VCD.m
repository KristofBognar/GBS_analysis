
% % load radnom VCD input file (all have sonde data up to 2012 from Cristen/Xiaoyi)
% load /home/kristof/work/GBS/PEARL-GBS/2017/VCD/vcd_input.mat

load /home/kristof/work/ozonesonde/Eureka/sonde_for_VCD.mat

clearvars -except sonde

yr=[];
doy=[];
ft=[];
o3=[];

start_yr=2018;

% load total columns from ozonesonde files
for year=start_yr:2019
    
    load(['/home/kristof/work/ozonesonde/Eureka/o3sonde_' num2str(year) '.mat'])
    
    % get year, day of year and fractional time
    formstr='yyyy-mm-dd HH:MM:SS';
    for i=1:length(f_list)
        datestr=[launchtime{i,1}, ' ', launchtime{i,2}];
        [ft_tmp,yr_tmp]=fracdate(datestr,formstr);
        
        % year is straightforward
        yr=[yr;yr_tmp];
        % day of year from fractional day
        doy=[doy; floor(ft_tmp)+1];
        % fractional time from fractional day
        ft=[ft; ft_tmp-floor(ft_tmp)];
        
        % save ozone column (in DU)
        o3=[o3;tot_col_DU(i)];
    end
    
end

% convert fractional time to hours
ft=ft*24;

% create one array
sonde2=[yr,doy,ft,o3];

% remove zeros from DU values (no corrected column from sonde)
sonde2(o3==0,:)=[];

% merge arrays
sonde(sonde(:,1)==2012,:)=[];

sonde=[sonde;sonde2];

clearvars -except sonde

save /home/kristof/work/ozonesonde/Eureka/sonde_for_VCD.mat


