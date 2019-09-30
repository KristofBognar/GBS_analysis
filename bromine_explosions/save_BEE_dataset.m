function save_BEE_dataset()
% Load and filter BrO profiles, and pair all available data
%
% high frequency observations are averaged for the time coverage of each profile
% low frequency observations are interpolated to the mean profile time


%% Load data

info=[];

% save some aerosol data
load('/home/kristof/work/profile_retrievals/profile_results/aerosol_profiles_filt_all.mat');

prof_aer=prof;
prof_aer_err=prof_err;
info_aer=info;

% BrO data
load('/home/kristof/work/profile_retrievals/profile_results/tracegas_profiles_filt_all.mat');

% BrO a priori info
load('/home/kristof/work/profile_retrievals/retr_times.mat');

% OPC data
load('/home/kristof/work/SMPS/smps+opc_2016-2018/opc_size_dist_all.mat');

% SMPS data
load('/home/kristof/work/SMPS/smps+opc_2016-2018/smps_size_dist_all.mat');

% PWS data
load('/home/kristof/work/weather_stations/ridge_lab/PWS_all.mat');
data((month(data.DateTime)>5 | month(data.DateTime)<3),:)=[];
pws_data=data;

% EWS data
load('/home/kristof/work/weather_stations/Eureka/EWS_PTU_and_weather_complete.mat')
ews_data=data;

% surface ozone data
load('/home/kristof/work/surface_ozone/surf_o3_all.mat');
surf_o3((month(surf_o3.DateTime)>5 | month(surf_o3.DateTime)<3),:)=[];

% pTOMCAT data
load('/home/kristof/work/models/pTOMCAT/pTOMCAT_ssa_all.mat');
ptom_time_ssa=ptom_time;
ptom_alt_ssa=ptom_alt;
load('/home/kristof/work/models/pTOMCAT/pTOMCAT_tg_all.mat');


%% Filter BrO profiles

ind_bad=find( info.DOFS<0.7 | info_aer.DOFS<0.7 | info_aer.col>5);
% percentiles with this filter (2015-2018):
% median: 9.1125e+12
% mean: 1.4578e+13
% 75: 1.5958e+13
% 80: 2.0782e+13
% 90: 3.4780e+13
% 95: 4.4484e+13

times(ind_bad)=[];
info(ind_bad,:)=[];
info_aer(ind_bad,:)=[];

prof(:,ind_bad)=[];
prof_err(:,ind_bad)=[];

prof_nd(:,ind_bad)=[];
prof_nd_err(:,ind_bad)=[];

prof_aer(:,ind_bad)=[];
prof_aer_err(:,ind_bad)=[];

avk(:,:,ind_bad)=[];

% get partial columns
part_prof=prof_nd*20000; % layers are 200m thick

% partial column height (top layer still included in the partial column)
top=4; % 600m layer, includes lab
% top=3; % 400m layer

below_lab=sum(part_prof(1:top,:))';
above_lab=sum(part_prof(top+1:end,:))';

% find daily time window used in profile retrieval

prof_len=match_prof_length(times);

prof_len=prof_len/2; % +- minutes around mean time


%% Pair all variables to the BrO profiles

%%% wind speed and direction
wspd_mean=find_coincident_mean(times, pws_data.DateTime, pws_data.WindSpd, prof_len);

wdir_mean=find_coincident_mean(times, pws_data.DateTime, pws_data.WindDir, prof_len, true);

%%% mean wind direction
% northerly winds, mean: 354 deg from gaussian fit, +-30 deg
ind_N=(~isnan(wspd_mean) & (wdir_mean >= 324 | wdir_mean < 24));
% southeasterly winds, mean: 123 deg from gaussian fit, +-30 deg  
ind_SE=(~isnan(wspd_mean) & (wdir_mean >= 93 & wdir_mean < 153));
% everything else
ind_rest=(~isnan(wspd_mean) & ( (wdir_mean < 93 & wdir_mean >= 24) |...
                              (wdir_mean < 324 & wdir_mean >= 153)) );

%%% PWS temperature, pressure, and pressure tendency from prev. hour
t_rl_mean=find_coincident_mean(times, pws_data.DateTime, pws_data.TempC, prof_len);
p_rl_mean=find_coincident_mean(times, pws_data.DateTime, pws_data.Pressmb, prof_len);
p_rl_mean_prev=find_coincident_mean(times-duration(1,0,0),...
                                    pws_data.DateTime, pws_data.Pressmb, 5);

%%% EWS temperature, pressure and winds
t_ews=interp1(ews_data.DateTime,ews_data.TempC,times,'linear','extrap');
p_ews=interp1(ews_data.DateTime,ews_data.StnPresskPa*10,times,'linear','extrap');
p_ews_prev=interp1(ews_data.DateTime,ews_data.StnPresskPa*10,times-duration(1,0,0),...
                   'linear','extrap');
wspd_ews=interp1(ews_data.DateTime,ews_data.WindSpdkmh*10/36,times,'linear','extrap');

tmp=unwrap(ews_data.WindDir10sdeg*10*pi/180)*180/pi;
tmp=interp1(ews_data.DateTime,tmp,times,'linear','extrap');
wdir_ews=mod(tmp,360);        

% cannot interpolate weather info column: find nearest indices
mjd2k_bro=mjd2k(times);
mjd2k_ews=mjd2k(ews_data.DateTime);
diff=abs(bsxfun(@minus,mjd2k_ews,mjd2k_bro'));
[~,ind]=min(diff,[],1);

ews_weather=ews_data.Weather(ind');

%%% T inversin from sonde
t_sonde=[];
wnd_sonde=[];
for yr=unique(times.Year)'
    t_sonde=[t_sonde; get_sonde_PT(yr)];
    wnd_sonde=[wnd_sonde; get_sonde_wnd(yr)];
end

t_0=interp1(t_sonde.date,t_sonde.T_0,times,'linear','extrap');
t_200=interp1(t_sonde.date,t_sonde.T_200,times,'linear','extrap');
t_400=interp1(t_sonde.date,t_sonde.T_400,times,'linear','extrap');
t_600=interp1(t_sonde.date,t_sonde.T_600,times,'linear','extrap');
t_inv_all=interp1(t_sonde.date,t_sonde.dT,times,'linear','extrap');
% t_inv_all=interp1(t_inv.date,t_inv.dT,times,'nearest','extrap');

wspd_0=interp1(wnd_sonde.date,wnd_sonde.wspd_0,times,'linear','extrap');
wspd_200=interp1(wnd_sonde.date,wnd_sonde.wspd_200,times,'linear','extrap');
wspd_400=interp1(wnd_sonde.date,wnd_sonde.wspd_400,times,'linear','extrap');
wspd_600=interp1(wnd_sonde.date,wnd_sonde.wspd_600,times,'linear','extrap');

tmp=unwrap(wnd_sonde.wdir_0*pi/180)*180/pi;
tmp=interp1(wnd_sonde.date,tmp,times,'linear','extrap');
wdir_0=mod(tmp,360);

tmp=unwrap(wnd_sonde.wdir_200*pi/180)*180/pi;
tmp=interp1(wnd_sonde.date,tmp,times,'linear','extrap');
wdir_200=mod(tmp,360);

tmp=unwrap(wnd_sonde.wdir_400*pi/180)*180/pi;
tmp=interp1(wnd_sonde.date,tmp,times,'linear','extrap');
wdir_400=mod(tmp,360);

tmp=unwrap(wnd_sonde.wdir_600*pi/180)*180/pi;
tmp=interp1(wnd_sonde.date,tmp,times,'linear','extrap');
wdir_600=mod(tmp,360);


%%% OPC data
% 1-10 microns (larger particles are extremely rare)
opc_sum=sum(opc_data(:,3:5),2);

ind=find(opc_sum==0 | isnan(opc_sum) | opc_sum>1.2); % there are a few spikes in the data
opc_data(ind,:)=[];
opc_time(ind)=[];
opc_tot_data(ind)=[];
opc_sum(ind)=[];

opc_mean=find_coincident_mean(times, opc_time, opc_sum, prof_len);

%%% SMPS data
% integrate Dp data
% 50-500 nm: 23-end
% 100-500 nm: 33-end
smps_sum=sum(smps_data(:,33:end),2);

ind=find(smps_sum==0 | isnan(smps_sum));
smps_data(ind,:)=[];
smps_time(ind)=[];
smps_tot_data(ind)=[];
smps_sum(ind)=[];

smps_mean=find_coincident_mean(times, smps_time, smps_sum, prof_len);

%%% surface ozone data
o3_mean=find_coincident_mean(times, surf_o3.DateTime, surf_o3.o3_ppb, prof_len);

%%% pTOMCAT data
% 0-4 km column
ptom_col_bro_out=interp1(ptom_time,ptom_col_bro,times,'linear','extrap');

% surface conc.
ptom_surf_bro_out=interp1(ptom_time,ptom_prof_bro(1,:),times,'linear','extrap');
ptom_surf_o3_out=interp1(ptom_time,ptom_prof_o3(1,:),times,'linear','extrap');

% supermicron aerosol profiles (both sea ice sourced and open ocean sourced)
% from ~1 micron to 10 micron
sism=sum(ptom_sissa(:,:,14:20),3); 
oosm=sum(ptom_oossa(:,:,14:20),3); 

% accumulation mode aerosol profiles (both sea ice sourced and open ocean sourced)
% from ~0.1 micron to 0.5 micron
sism_accum=sum(ptom_sissa(:,:,7:11),3); 
oosm_accum=sum(ptom_oossa(:,:,7:11),3); 

% interpolate aerosol to ridge lab altitude
sism_rl=NaN(1,length(ptom_time_ssa));
oosm_rl=NaN(1,length(ptom_time_ssa));
sism_rl_accum=NaN(1,length(ptom_time_ssa));
oosm_rl_accum=NaN(1,length(ptom_time_ssa));

for i=1:length(ptom_time_ssa)
    sism_rl(i)=interp1(ptom_alt_ssa(:,i),sism(:,i),610,'linear');
    oosm_rl(i)=interp1(ptom_alt_ssa(:,i),oosm(:,i),610,'linear');
    sism_rl_accum(i)=interp1(ptom_alt_ssa(:,i),sism_accum(:,i),610,'linear');
    oosm_rl_accum(i)=interp1(ptom_alt_ssa(:,i),oosm_accum(:,i),610,'linear');
end

% interpolate aerosol to BrO profile times
ptom_sissa=interp1(ptom_time_ssa,sism_rl,times,'linear','extrap');
ptom_oossa=interp1(ptom_time_ssa,oosm_rl,times,'linear','extrap');
ptom_sissa_accum=interp1(ptom_time_ssa,sism_rl_accum,times,'linear','extrap');
ptom_oossa_accum=interp1(ptom_time_ssa,oosm_rl_accum,times,'linear','extrap');

% T and P from pTOMCAT data
% interpolate to ridge lab altitude
ptom_T_0_tmp=ptom_T(1,:);
ptom_T_200_tmp=NaN(1,length(ptom_time));
ptom_T_600_tmp=NaN(1,length(ptom_time));

ptom_P_0_tmp=ptom_P(1,:);
ptom_P_600_tmp=NaN(1,length(ptom_time));

for i=1:length(ptom_time)
    ptom_T_200_tmp(i)=interp1(ptom_alt(:,i),ptom_T(:,i),610,'linear');
    ptom_T_600_tmp(i)=interp1(ptom_alt(:,i),ptom_T(:,i),610,'linear');
    ptom_P_600_tmp(i)=interp1(ptom_alt(:,i),ptom_P(:,i),610,'linear');
end

% interpolate t and P to BrO profile times
ptom_T_0=interp1(ptom_time,ptom_T_0_tmp,times,'linear','extrap')-273.15;
ptom_T_200=interp1(ptom_time,ptom_T_200_tmp,times,'linear','extrap')-273.15;
ptom_T_600=interp1(ptom_time,ptom_T_600_tmp,times,'linear','extrap')-273.15;

ptom_P_0=interp1(ptom_time,ptom_P_0_tmp,times,'linear','extrap')*1e-2;
ptom_P_0_prev=interp1(ptom_time,ptom_P_0_tmp,times-duration(1,0,0),...
                      'linear','extrap')*1e-2;
ptom_P_600=interp1(ptom_time,ptom_P_600_tmp,times,'linear','extrap')*1e-2;
ptom_P_600_prev=interp1(ptom_time,ptom_P_600_tmp,times-duration(1,0,0),...
                        'linear','extrap')*1e-2;


%% Create and save table

bee_dataset=table();

bee_dataset.times=times;

bee_dataset.bro_col=info.col;
bee_dataset.bro_col_err=info.col_err;
bee_dataset.bro_surf_ppt=prof(1,:)'*1e6;
bee_dataset.bro_surf_ppt_err=prof_err(1,:)'*1e6;
bee_dataset.bro_col_ratio=below_lab./info.col;

bee_dataset.aer_ext=info_aer.col;

bee_dataset.wspd_ms=wspd_mean;
bee_dataset.wdir=wdir_mean;

bee_dataset.N_SE_rest=NaN(size(bee_dataset.bro_col));
bee_dataset.N_SE_rest(ind_N)=1;
bee_dataset.N_SE_rest(ind_SE)=2;
bee_dataset.N_SE_rest(ind_rest)=3;

bee_dataset.wspd_ms_EWS=wspd_ews;
bee_dataset.wdir_EWS=wdir_ews;
bee_dataset.EWS_obs=ews_weather;

bee_dataset.T_PWS=t_rl_mean;
bee_dataset.T_EWS=t_ews;

bee_dataset.P_PWS=p_rl_mean;
bee_dataset.P_PWS_tend=p_rl_mean-p_rl_mean_prev;
bee_dataset.P_EWS=p_ews;
bee_dataset.P_EWS_tend=p_ews-p_ews_prev;

bee_dataset.sonde_T_0=t_0;
bee_dataset.sonde_T_200=t_200;
bee_dataset.sonde_T_400=t_400;
bee_dataset.sonde_T_600=t_600;
bee_dataset.sonde_dT=t_inv_all;

bee_dataset.sonde_wspd_0=wspd_0;
bee_dataset.sonde_wspd_200=wspd_200;
bee_dataset.sonde_wspd_400=wspd_400;
bee_dataset.sonde_wspd_600=wspd_600;

bee_dataset.sonde_wdir_0=wdir_0;
bee_dataset.sonde_wdir_200=wdir_200;
bee_dataset.sonde_wdir_400=wdir_400;
bee_dataset.sonde_wdir_600=wdir_600;

bee_dataset.OPC_supermicron=opc_mean;

bee_dataset.SMPS_100_500=smps_mean;

bee_dataset.o3_surf=o3_mean;

bee_dataset.bro_col_ptom=ptom_col_bro_out;
bee_dataset.bro_surf_ptom=ptom_surf_bro_out;
bee_dataset.o3_surf_ptom=ptom_surf_o3_out;
bee_dataset.ptom_supermicron=ptom_sissa+ptom_oossa;
bee_dataset.ptom_supermicron_si_frac=ptom_sissa./(ptom_sissa+ptom_oossa);
bee_dataset.ptom_submicron=ptom_sissa_accum+ptom_oossa_accum;
bee_dataset.ptom_submicron_si_frac=ptom_sissa_accum./(ptom_sissa_accum+ptom_oossa_accum);

bee_dataset.ptom_T_0=ptom_T_0;
bee_dataset.ptom_T_200=ptom_T_200;
bee_dataset.ptom_T_600=ptom_T_600;
bee_dataset.ptom_P_0=ptom_P_0;
bee_dataset.ptom_P_0_tend=ptom_P_0-ptom_P_0_prev;
bee_dataset.ptom_P_600=ptom_P_600;
bee_dataset.ptom_P_600_tend=ptom_P_600-ptom_P_600_prev;


save('/home/kristof/work/BEEs/BEE_dataset_all.mat','bee_dataset')



end

function prof_len=match_prof_length(times)

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




