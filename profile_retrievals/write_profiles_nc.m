% function write_profiles_nc(option, year)
% write aerosol or BrO profiles to netCDF file

year=2019;
option='tg';

% directories and output file name
if option=='tg';
    
    species_str='BrO';
    
    savedir=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
             num2str(year) '/'];
%     savedir=['/home/kristof/work/profile_retrievals/profile_results/'];

    f_in=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
          num2str(year) '/tracegas/' 'profiles_' num2str(year) '_filt.mat'];

    f_out=[savedir 'Eureka_BrO_profiles_' num2str(year) '.nc'];

else
    
    species_str='aerosol';
    
    savedir=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
             num2str(year) '/'];

    f_in=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
          num2str(year) '/aerosol/' 'profiles_' num2str(year) '_filt.mat'];

    f_out=[savedir 'Eureka_aerorol_profiles_' num2str(year) '.nc'];
    
end

% load profiles
% f_in=[savedir 'tracegas_profiles_' num2str(year) '_normal_aer_only.mat'];

load(f_in)

if exist(f_out,'file'), delete(f_out); end
ncid = netcdf.create(f_out,'netcdf4'); % change type to netcdf4 (default is netcdf4_classic)
netcdf.close(ncid);

% create netCDF file with variables
nccreate(f_out, 'altitude','Dimensions',{'Altitude',length(alt)});
nccreate(f_out, 'time','Dimensions',{'Time',length(ft)});

% write tracegas-specific data
if option=='tg'
    nccreate(f_out, 'profile_vmr','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
    nccreate(f_out, 'profile_vmr_error','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
    nccreate(f_out, 'profile_nd','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
    nccreate(f_out, 'profile_nd_error','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
else
    nccreate(f_out, 'profile','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
    nccreate(f_out, 'profile_error','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
end
nccreate(f_out, 'column','Dimensions',{'Time',length(ft)});
nccreate(f_out, 'column_error','Dimensions',{'Time',length(ft)});
nccreate(f_out, 'dofs','Dimensions',{'Time',length(ft)});
nccreate(f_out, 'avk_col','Dimensions',{'Altitude',length(alt),'Time',length(ft)});
nccreate(f_out, 'avk','Dimensions',{'Altitude',length(alt),'Altitude',length(alt),'Time',length(ft)});


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

% Add variable attributes
ncwriteatt(f_out,'/altitude','units','km');
ncwriteatt(f_out,'/time','units','mjd2k, zero at Jan. 1, 2000 at 00:00 (UTC)');

if option=='tg'
    ncwriteatt(f_out,'/profile_vmr','units','parts per trillion (ppt)');
    ncwriteatt(f_out,'/profile_vmr_error','units','parts per trillion (ppt)');
    ncwriteatt(f_out,'/profile_nd','units','molecules per cm^3');
    ncwriteatt(f_out,'/profile_nd_error','units','molecules per cm^3');
    ncwriteatt(f_out,'/column','units','molecules per cm^2');
    ncwriteatt(f_out,'/column_error','units','molecules per cm^2');
else
    ncwriteatt(f_out,'/profile','units','extinction, km^-1');
    ncwriteatt(f_out,'/profile_error','units','extinction, km^-1');
    ncwriteatt(f_out,'/column','units','aerosol optical depth');
    ncwriteatt(f_out,'/column_error','units','aerosol optical depth');
end    
ncwriteatt(f_out,'/dofs','description','degrees of freedom for signal for each profile');
ncwriteatt(f_out,'/avk_col','description','column averaging kernel');
ncwriteatt(f_out,'/avk','description','full averaging kernel, altitude x altitude specific avk x time');

% write variable values
ncwrite(f_out, 'time',ft_to_mjd2k(ft,year));
ncwrite(f_out, 'altitude',alt);
if option=='tg'
    ncwrite(f_out, 'profile_vmr',prof*1e6);
    ncwrite(f_out, 'profile_vmr_error',prof_err*1e6);
    ncwrite(f_out, 'profile_nd', prof_nd*1e6);
    ncwrite(f_out, 'profile_nd_error',prof_nd_err*1e6);
else
    ncwrite(f_out, 'profile',prof);
    ncwrite(f_out, 'profile_error',prof_err);
end
ncwrite(f_out, 'column',info.col);
ncwrite(f_out, 'column_error',info.col_err);
ncwrite(f_out, 'dofs',info.DOFS);
ncwrite(f_out, 'avk_col',avk_col);
ncwrite(f_out, 'avk',avk);


% end


%% write BrO profiles as text file

if option=='tg'
    
    prof=prof';
    prof_err=prof_err';
    prof_nd=prof_nd';
    prof_nd_err=prof_nd_err';

    for i=1:length(alt)

        prof_head{i}=['prof_' num2str(alt(i)*1000) 'm'];
        prof_err_head{i}=['prof_err_' num2str(alt(i)*1000) 'm'];

    end

    out=array2table([ft,info.DOFS,info.col,info.col_err,prof*1e6,prof_err*1e6],...
                    'VariableNames',[{'Fractional_time','DOFS','column','column_err'},...
                    prof_head, prof_err_head]);

    writetable(out,['Eureka_BrO_profiles_ppt_' num2str(year)],'delimiter',',');

% %     add manually to file:
% %     # Fractional time is 0 at Jan 1, 00:00
% %     # Columns are 0-4 km, in molec/cm^2
% %     # profiles and errors are parts per trillion by volume

end
