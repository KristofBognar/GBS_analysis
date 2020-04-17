function write_profiles_nc(option, times, alt, prof, prof_err, info, avk, avk_col,...
                           prof_nd, prof_nd_err)
% write aerosol or BrO profiles to netCDF file

% get time info (still fractional time, MJD2k conversion happens later)
[ft,year]=fracdate(times);

% check if more than one year
if length(unique(year))==1
    year=year(1);
%     ft=ft_to_mjd2k(ft,year);
else
    error('Input must contain one year only')
end

% set output directory
savedir=['/home/kristof/work/profile_retrievals/profile_results/nc_files/'];
    
% directories and output file name
if strcmp(option,'tg');
    
    species_str='BrO';
%     f_in=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
%           num2str(year) '/tracegas/' 'profiles_' num2str(year) '_filt.mat'];

elseif strcmp(option,'aer');
    
    species_str='aerosol';
%     f_in=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
%           num2str(year) '/aerosol/' 'profiles_' num2str(year) '_filt.mat'];
    
end

f_out=[savedir 'Eureka_MAX-DOAS_' species_str '_profiles_' ...
                datestr(ft_to_date(ft(1),year),'yyyymmddTHHMMSSZ') '_' ...
                datestr(ft_to_date(ft(end),year),'yyyymmddTHHMMSSZ') '.nc'];

% load profiles
% f_in=[savedir 'tracegas_profiles_' num2str(year) '_normal_aer_only.mat'];
% load(f_in)

%% create file
if exist(f_out,'file'), delete(f_out); end
ncid = netcdf.create(f_out,'netcdf4'); % change type to netcdf4 (default is netcdf4_classic)
netcdf.close(ncid);


%% create variables
% create netCDF file with variables
nccreate(f_out, 'altitude','Dimensions',{'Altitude',length(alt)});
nccreate(f_out, 'time','Dimensions',{'Time',length(ft)});

% write tracegas-specific data
if strcmp(option,'tg');
    nccreate(f_out, 'profile_vmr','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
    nccreate(f_out, 'profile_vmr_error','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
    nccreate(f_out, 'profile_nd','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
    nccreate(f_out, 'profile_nd_error','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
elseif strcmp(option,'aer');
    nccreate(f_out, 'profile','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
    nccreate(f_out, 'profile_error','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
end

nccreate(f_out, 'column','Dimensions',{'Time',length(ft)});
nccreate(f_out, 'column_error','Dimensions',{'Time',length(ft)});
nccreate(f_out, 'dofs','Dimensions',{'Time',length(ft)});
nccreate(f_out, 'avk_col','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
nccreate(f_out, 'avk','Dimensions',{'Altitude',length(alt),'Altitude',length(alt),'Time',length(ft)});

nccreate(f_out, 'scan_length','Dimensions',{'Time',length(ft)});
nccreate(f_out, 'apriori_surf','Dimensions',{'Time',length(ft)});
nccreate(f_out, 'apriori_h','Dimensions',{'Time',length(ft)});


% add global attributes
ncwriteatt(f_out,'/','summary',['MAX-DOAS profiles of '...
                                species_str ' using the HEIPro retrievl algorithm']);
ncwriteatt(f_out,'/','instrument','PEARL-GBS');
% ncwriteatt(f_out,'/','elevation_angles','90,30,15,10,5,2,1,-1');
ncwriteatt(f_out,'/','retrieval_wavelength','360.8 nm');
ncwriteatt(f_out,'/','location','PEARL Ridge Lab, Eureka, Nunavut, Canada');
ncwriteatt(f_out,'/','location_lat','80.053');
ncwriteatt(f_out,'/','location_lon','-86.416');
ncwriteatt(f_out,'/','location_alt','0.610');
ncwriteatt(f_out,'/','start_date',datestr(ft_to_date(ft(1),year),'yyyymmddTHHMMSSZ'));
ncwriteatt(f_out,'/','end_date',datestr(ft_to_date(ft(end),year),'yyyymmddTHHMMSSZ'));
ncwriteatt(f_out,'/','creation_date',datestr(datetime('now','timezone','utc'),...
                                             'yyyymmddTHHMMSSZ'));
ncwriteatt(f_out,'/','creator_name','Kristof Bognar');
ncwriteatt(f_out,'/','institution','Department of Physics, University of Toronto');
ncwriteatt(f_out,'/','project_PI','Kimberly Strong');

% add variable attributes
ncwriteatt(f_out,'/altitude','units','km');
ncwriteatt(f_out,'/time','units','mjd2k, days since 00:00, Jan. 1, 2000 (UTC) (12:00, Jan. 1, 2000 = 0.5)');
ncwriteatt(f_out,'/time','description','mean time of each profile');

if strcmp(option,'tg');
    ncwriteatt(f_out,'/profile_vmr','units','parts per trillion (ppt)');
    ncwriteatt(f_out,'/profile_vmr_error','units','parts per trillion (ppt)');
    ncwriteatt(f_out,'/profile_nd','units','molecules per cm^3');
    ncwriteatt(f_out,'/profile_nd_error','units','molecules per cm^3');
    ncwriteatt(f_out,'/column','units','molecules per cm^2');
    ncwriteatt(f_out,'/column_error','units','molecules per cm^2');
    
    ncwriteatt(f_out,'/apriori_surf','units','parts per trillion (ppt)');

elseif strcmp(option,'aer');
    ncwriteatt(f_out,'/profile','units','extinction, km^-1');
    ncwriteatt(f_out,'/profile_error','units','extinction, km^-1');
    ncwriteatt(f_out,'/column','units','aerosol optical depth');
    ncwriteatt(f_out,'/column_error','units','aerosol optical depth');
    
    ncwriteatt(f_out,'/apriori_surf','units','extinction, km^-1');
    
end 

ncwriteatt(f_out,'/dofs','description','degrees of freedom for signal for each profile');
ncwriteatt(f_out,'/avk_col','description','column averaging kernel');
ncwriteatt(f_out,'/avk','description','full averaging kernel, altitude x altitude specific avk x time');

ncwriteatt(f_out,'/apriori_surf','description','surface value of the exponential a priori profile');
ncwriteatt(f_out,'/apriori_h','units','km');
ncwriteatt(f_out,'/apriori_h','description','scale height of the exponential a priori profile');
ncwriteatt(f_out,'/scan_length','units','minutes');
ncwriteatt(f_out,'/scan_length','description','total duration of each profile (duration of the corresponding MAX-DOAS scan)');



%% write variable values
ncwrite(f_out, 'time',ft_to_mjd2k(ft,year));
ncwrite(f_out, 'altitude',alt);
if strcmp(option,'tg');
    ncwrite(f_out, 'profile_vmr',prof*1e6);
    ncwrite(f_out, 'profile_vmr_error',prof_err*1e6);
    ncwrite(f_out, 'profile_nd', prof_nd*1e6);
    ncwrite(f_out, 'profile_nd_error',prof_nd_err*1e6);
    
    
elseif strcmp(option,'aer');
    ncwrite(f_out, 'profile',prof);
    ncwrite(f_out, 'profile_error',prof_err);
end
ncwrite(f_out, 'column',info.col);
ncwrite(f_out, 'column_error',info.col_err);
ncwrite(f_out, 'dofs',info.DOFS);
ncwrite(f_out, 'avk_col',avk_col);
ncwrite(f_out, 'avk',avk);

ncwrite(f_out, 'scan_length',info.prof_len);
ncwrite(f_out, 'apriori_surf',info.ap_surf);
ncwrite(f_out, 'apriori_h',info.ap_h);

%% check for aux data input   

if ismember('temperature', info.Properties.VariableNames)
    nccreate(f_out, 'temperature','Dimensions',{'Time',length(ft)});
    ncwriteatt(f_out,'/temperature','units','K');
    ncwriteatt(f_out,'/temperature','description','mean temperature for the profile duration (measured at the PEARL Ridge Lab)');
    ncwriteatt(f_out,'/temperature','fill_value','-9999');
    ncwrite(f_out, 'temperature',info.temperature);
end

if ismember('pressure', info.Properties.VariableNames)
    nccreate(f_out, 'pressure','Dimensions',{'Time',length(ft)});
    ncwriteatt(f_out,'/pressure','units','Pa');
    ncwriteatt(f_out,'/pressure','description','mean pressure for the profile duration (measured at the PEARL Ridge Lab)');
    ncwriteatt(f_out,'/pressure','fill_value','-9999');
    ncwrite(f_out, 'pressure',info.pressure);
end

if ismember('wind_speed', info.Properties.VariableNames)
    nccreate(f_out, 'wind_speed','Dimensions',{'Time',length(ft)});
    ncwriteatt(f_out,'/wind_speed','units','m s^-1');
    ncwriteatt(f_out,'/wind_speed','description','mean wind speed for the profile duration (measured at the PEARL Ridge Lab)');
    ncwriteatt(f_out,'/wind_speed','fill_value','-9999');
    ncwrite(f_out, 'wind_speed',info.wind_speed);
end

if ismember('wind_dir', info.Properties.VariableNames)
    nccreate(f_out, 'wind_dir','Dimensions',{'Time',length(ft)});
    ncwriteatt(f_out,'/wind_dir','units','degrees, N=0, E=90');
    ncwriteatt(f_out,'/wind_dir','description','mean wind direction for the profile duration (measured at the PEARL Ridge Lab)');
    ncwriteatt(f_out,'/wind_dir','fill_value','-9999');
    ncwrite(f_out, 'wind_dir',info.wind_dir);
end
    

end


%% write BrO profiles as text file

% if option=='tg'
%     
%     prof=prof';
%     prof_err=prof_err';
%     prof_nd=prof_nd';
%     prof_nd_err=prof_nd_err';
% 
%     for i=1:length(alt)
% 
%         prof_head{i}=['prof_' num2str(alt(i)*1000) 'm'];
%         prof_err_head{i}=['prof_err_' num2str(alt(i)*1000) 'm'];
% 
%     end
% 
%     out=array2table([ft,info.DOFS,info.col,info.col_err,prof*1e6,prof_err*1e6],...
%                     'VariableNames',[{'Fractional_time','DOFS','column','column_err'},...
%                     prof_head, prof_err_head]);
% 
%     writetable(out,['Eureka_BrO_profiles_ppt_' num2str(year)],'delimiter',',');
% 
% % %     add manually to file:
% % %     # Fractional time is 0 at Jan 1, 00:00
% % %     # Columns are 0-4 km, in molec/cm^2
% % %     # profiles and errors are parts per trillion by volume
% 
% end
