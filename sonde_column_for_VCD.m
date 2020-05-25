function sonde_column_for_VCD(start_yr, end_year)
%%% function to generate ozone column inputs required for VCD retrieval
%%% code loads ozonesonde data saved by read_ozonesonde.m
%%%

% load old VCD input file 
if ismac
    error('Get sample VCD file and set up path')
elseif isunix

    file_path='/home/kristof/work/ozonesonde/Eureka/';
    sonde_path='/home/kristof/work/ozonesonde/Eureka/';
    
    load([file_path 'sonde_for_VCD.mat'])
end

yr=[];
doy=[];
ft=[];
o3=[];

% start_yr=2018;
% end_year=2020;

% load total columns from ozonesonde files
for year=start_yr:end_year
    
    load([sonde_path 'o3sonde_' num2str(year) '.mat'])
    
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
    
    % check for overlaping year in sonde file
    % don't overwrite anything, ask user to handle issue
    if sum(sonde(:,1)==year)>0,
        error(['Data from ' num2str(year) ' already present in sonde file:\n%s'],...
              'Adjust start year or delete partial year from the saved file')
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

save([file_path 'sonde_for_VCD.mat'],'sonde')

end
