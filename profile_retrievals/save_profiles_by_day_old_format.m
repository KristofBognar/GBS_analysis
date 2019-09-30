
% plot and save retrieved aerosol and tracegas profiles

%% input parameters

year=2018;

% profiles to plot
% 'a' for aerosol, 'tg' for tracegas
option='a'; 

% select directories and profile column numbers
if option=='a'

    ymd='20180329';
    version='noadaptive_skip0';
    dir_name=[ymd '_' version];
    if any(year==[2017,2018]), dir_name=ymd; end

    prof_dir=['/home/kristof/Drive/Retrieval_settings_A (3)/',dir_name,'/general/'];
    savedir=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
             num2str(year) '/aerosol/'];

elseif option=='tg'

%     units='vmr'; % 'vmr', 'cnc', or 'nd'
                 % vmr in ppm, cnc in mu g /m^3, nd in moled/m^3
    
    ymd='20160322';
    version='noadaptive_1iter_skip0';
    dir_name=[ymd '_' version];
    if any(year==[2017,2018]), dir_name=ymd; end

    prof_dir=['/home/kristof/Drive/Retrieval_settings_A (2) (2)/',dir_name,'/general/'];
    savedir=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
             num2str(year) '/tracegas/'];

else
    error('select aerosol or tracegas files to read')
end

cur_dir=(pwd);

%% read filenames and get profile times

cd(prof_dir)

if option=='tg'
    %% read number density profiles
    units='nd';
    % profile file name
    search=['*' units '_prof_20*' ];
    tmp = dir(search); 
    f_prof=tmp.name;
    % profile error file
    search=['*' units '_prof_errs*' ];
    tmp = dir(search); 
    f_prof_err=tmp.name;
    
    % read profiles and errors
    prof_tmp=dlmread(f_prof,'',1,0);
    alt=prof_tmp(:,1);
    prof_nd=prof_tmp(:,2:end);

    prof_err_tmp=dlmread(f_prof_err,'',1,0);
    prof_nd_err=prof_err_tmp(:,2:end);

    %% read vmr profiles
    units='vmr';
    % profile file name
    search=['*' units '_prof_20*' ];
    tmp = dir(search); 
    f_prof=tmp.name;
    % profile error file
    search=['*' units '_prof_errs*' ];
    tmp = dir(search); 
    f_prof_err=tmp.name;
    
    % read profiles and errors
    prof_tmp=dlmread(f_prof,'',1,0);
    alt=prof_tmp(:,1);
    prof=prof_tmp(:,2:end);

    prof_err_tmp=dlmread(f_prof_err,'',1,0);
    prof_err=prof_err_tmp(:,2:end);

    
elseif option=='a'
    % profile file name
    tmp = dir('all_profiles_*'); 
    f_prof=tmp.name;
    % profile error file
    tmp = dir('all_prof_errs_*'); 
    f_prof_err=tmp.name;
  
    % read profiles and errors
    prof_tmp=dlmread(f_prof,'',1,0);
    alt=prof_tmp(:,1);
    prof=prof_tmp(:,2:end);

    prof_err_tmp=dlmread(f_prof_err,'',1,0);
    prof_err=prof_err_tmp(:,2:end);
    
    % dummy nd profile variables so I can use the script for aerosols as well
    prof_nd=0;
    prof_nd_err=0;
end    

%% read measured/retrieved DSCDs
tmp = dir('meas_*'); 
f_dscd=tmp.name;

dscd=dlmread(f_dscd,'',1,2);

%% read profile header and get ft

% read header
fid=fopen(f_prof,'r');
tmp=textscan(fid,'%s',1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
fclose(fid);

% get time info
tmp=strsplit(tmp{1}{1},' ');

times=tmp(3:end);

for i=1:size(times,2)
    
    times{i}=[ymd, '_', times{i}];
end
    
formstr='yyyymmdd_HH:MM:SS';
[ft]=fracdate(times, formstr);

%% kick out bad profiles

ind=find(prof(1,:)==-9999);
if ~isempty(ind), ft(ind)=[]; end

prof(prof==-9999)=NaN;
prof_err(prof_err==-9999)=NaN;

prof=prof(:,all(~isnan(prof)));
prof_err=prof_err(:,all(~isnan(prof_err)));

% same for nd profile
% columns that failed should be the same as for vmr profile
prof_nd(prof_nd==-9999)=NaN;
prof_nd_err(prof_nd_err==-9999)=NaN;

prof_nd=prof_nd(:,all(~isnan(prof_nd)));
prof_nd_err=prof_nd_err(:,all(~isnan(prof_nd_err)));

%% read DSCD time columns and get ft
fid=fopen(f_dscd,'r');
tmp = textscan(fid,'%10s%9s%[^\n\r]', 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fid);

% Remove white space around all cell columns.
tmp{1} = strtrim(tmp{1});
tmp{2} = strtrim(tmp{2});
times_dscd=cell(size(tmp{1}));

for i=1:size(tmp{1},1)
    times_dscd{i}=[tmp{1}{i} ' ' tmp{2}{i}];
end

formstr='dd/mm/yyyy HH:MM:SS';
[ft_dscd]=fracdate(times_dscd, formstr);

%% read retrieval details

cd('../retrieval_details')

% get filenames
tmp = dir('retr_*'); 
f_info={tmp.name}';

% if max(size(f_info))~=size(ft,1), error('time mismatch'), end

% define arrays to hold values
info=zeros(size(ft,1),4); 

for i=1:max(size(f_info))

    % read all data
    tmp=dlmread(f_info{i},':',0,1);

    % save DoF
    info(i,1)=tmp(4);
    % save VCD (molec/cm^2) for tg
    % save optical thickness for aer
    info(i,2)=tmp(7);
    % save VCD error
    % save optical thickness error
    info(i,3)=tmp(8);

    % get time info
    tmp=strsplit(f_info{i},'.');
    % split along underscores
    tmp=strsplit(tmp{1},'_');
    % get date and time
    times_info=[tmp{end-1},'_',tmp{end}];
    
    formstr='yyyymmdd_HHMM';
    [ft_info]=fracdate(times_info, formstr);
    
    info(i,4)=ft_info;

end

cd(cur_dir)

%% plot profiles

figure(1)

% color plot of profiles during the day
surf(ft,alt,prof,'EdgeColor','None', 'facecolor', 'interp')
view(2)
colormap(jet(300))
colorbar
% xlim([ft(1),ft(end)])
ylim([-0.1,4.1])

xlabel(['Fractional day, ' num2str(year) ' (UTC)'])
ylabel('Altitude (km)')

% % %% plot measured/modelled DSCDs
% % 
% % figure(2)
% % 
% % plot(ft_dscd,dscd(:,5),'bo'), hold on
% % plot(ft_dscd,dscd(:,end),'rx')

%% save results

f_out=[savedir dir_name '.mat'];

if ~exist(savedir, 'dir'), mkdir(savedir); end

clearvars -except alt prof prof_nd prof_err prof_nd_err ft dscd ft_dscd f_out info

save(f_out)
