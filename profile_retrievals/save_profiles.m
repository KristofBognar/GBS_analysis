
% plot and save retrieved aerosol and tracegas profiles
% This code saves results by day, use Merge_profiles.m to create yearly files

%% input parameters

% profiles to read
% 'a' for aerosol, 'tg' for tracegas
option='tg'; 

% year=2018;
year=1; % for reading tests/rerurns for multiple years

% elevation correction used in the retrieval (to correctly read in number of elevs used)
elev_corr=0;

% location, 'Eureka' or 'CINDI-2'
loc='eureka';

% string after ymd in folder names, if applicable
folder_str='';

% default results directory, optput folders renamed just eureka_<year>
% other versions: output folders named  just eureka_<year>_<version>
default_version_eureka='Retrieval_settings_A';

% version_eureka=default_version_eureka;
% version_eureka='aer_10iter';
% version_eureka='surf_ppt_10'; % for reading a priori test results
% version_eureka='surf_ppt_5to1'; % surf conc = 5 ppt replaced with 1 ppt
% version_eureka='surf_ppt_1to5'; % surf conc = 1 ppt replaced with 5 ppt
version_eureka='det_lim_test'; % august 2018, for detection limit



%% years and dates
if year==2010
    % only retrieved select dates with promising BEEs
    start_date=datetime(year,4,10,12,0,0);
    end_date=datetime(year,5,20,12,0,0);
    
elseif year==2011
    % % only retrieved dates used in Zhao et al., 2016
    start_date=datetime(year,4,1,12,0,0);
    end_date=datetime(year,4,5,12,0,0);
    elev_corr=1;
    
elseif year==2013
    % april 3 missing, March 19, 28-30, April 2,4,6,7,9 have crazy aer values,
    % so BrO retrieval fails (dates removed from aer and BrO results)
    start_date=datetime(year,3,10,12,0,0);
    end_date=datetime(year,4,23,12,0,0);
    
elseif year==2015
    start_date=datetime(year,3,3,12,0,0);
    end_date=datetime(year,5,31,12,0,0);
    
elseif year==2016
    % datamissing in March due to cold weekend, and tracker misalignment (past 23rd)
    % shipped shouth for CINDI-2 in May
    start_date=datetime(year,3,6,12,0,0);
    end_date=datetime(year,5,11,12,0,0);
    
elseif year==2017
    % tracker stoped working in April
    start_date=datetime(year,3,7,12,0,0);
    end_date=datetime(year,4,04,12,0,0);
    
elseif year==2018
    % May 13 is missing
    start_date=datetime(year,3,5,12,0,0);
    end_date=datetime(year,5,31,12,0,0);
%     start_date=datetime(year,3,20,12,0,0);
%     end_date=datetime(year,4,20,12,0,0);

elseif year==2019
    % April 22, 28 missing
    start_date=datetime(year,3,5,12,0,0);
    end_date=datetime(year,5,31,12,0,0);
    
elseif year==1
    % reading a priori test results
    start_date=datetime(2015,3,1,12,0,0);
    end_date=datetime(2019,5,31,12,0,0);
    
end

% cindi-2 retrieval tests
% year=2016;
% start_date=datetime(year,9,12,12,0,0);
% end_date=datetime(year,9,28,12,0,0);
% % folder_str='_aer_20s'; % _fW,_20s,_fW_20s,_strat_aer,_udo_all,_udo_ap, + aer_*,_1i_20s for no2
% loc='CINDI-2';
% version_cindi2='medianDSCD_vis_v3_all';



% plot individual days?
show_plot=0;

% give warning if elevation correction is set
if elev_corr
   
    disp(['Warning: elevation correction set to ' num2str(elev_corr) ' deg. Proceed?'])
    disp('(return): continue; (n): reset elev corr to 0')
    tmp=input('','s');
    
    if ~isempty(tmp)
        elev_corr=0;
        disp('Elev corr reset to 0')
    end
    
end

%% loop over all days
for current_date=start_date:1:end_date

    % current day string (folder names)
    ymd=datestr(current_date,'yyyymmdd');
    
    %% select directories
    if option=='a'

        if strcmp(loc,'eureka')
            prof_dir=['/media/kristof/Windows7_OS/SCIATRAN2/AEROSOL_RETRIEVAL_v-1-2/',...
                       'Campaign/' version_eureka '/',ymd,folder_str,'/general/'];

            if strcmp(version_eureka,default_version_eureka)
                savedir=['/home/kristof/work/profile_retrievals/profile_results/',loc,'_'...
                         num2str(year) '/aerosol/'];
            else
                savedir=['/home/kristof/work/profile_retrievals/profile_results/',loc,'_'...
                         num2str(year) '_' version_eureka '/aerosol/'];
            end 
            
        elseif strcmp(loc,'CINDI-2')         
            prof_dir=['/media/kristof/Windows7_OS/SCIATRAN2/AEROSOL_RETRIEVAL_v-1-2/',...
                       'Campaign/Retrieval_settings_A/',version_cindi2,'/',ymd,folder_str,...
                       '/general/'];

            savedir=['/home/kristof/work/profile_retrievals/profile_results/',...
                      'CINDI-2/matlab_files/',version_cindi2,'/aerosol/'];
        end 

    elseif option=='tg'

        if strcmp(loc,'eureka')
            prof_dir=['/media/kristof/Windows7_OS/SCIATRAN2/TRACEGAS_RETRIEVAL_v-1-2/',...
                       'Campaign/' version_eureka '/',ymd,folder_str,'/general/'];

            if strcmp(version_eureka,default_version_eureka)
                savedir=['/home/kristof/work/profile_retrievals/profile_results/',loc,'_'...
                         num2str(year) '/tracegas/'];
            else
                savedir=['/home/kristof/work/profile_retrievals/profile_results/',loc,'_'...
                         num2str(year) '_' version_eureka '/tracegas/'];
            end
            
        elseif strcmp(loc,'CINDI-2')
            prof_dir=['/media/kristof/Windows7_OS/SCIATRAN2/TRACEGAS_RETRIEVAL_v-1-2/',...
                       'Campaign/Retrieval_settings_A/',version_cindi2,'/',ymd,folder_str,...
                       '/general/'];

            savedir=['/home/kristof/work/profile_retrievals/profile_results/',...
                      'CINDI-2/matlab_files/',version_cindi2,'/tracegas/'];
        end
    end

    cur_dir=(pwd);

    %% read filenames and get profile times

    % skip mising days
    try
        cd(prof_dir)
    catch
        continue
    end
    
    if option=='tg'
        %% read number density profiles
        units='nd'; % vmr in ppm, cnc in mu g /m^3, nd in moled/m^3
        
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
        units='vmr';% vmr in ppm, cnc in mu g /m^3, nd in moled/m^3
        
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
        prof_nd=prof;
        prof_nd_err=prof_err;
        
    end    

    %% read measured/retrieved DSCDs
    tmp = dir('meas_*'); 
    f_dscd=tmp.name;

    dscd=read_prof_dscd(f_dscd);
    
    % save dscd times as fractional time
    ft_dscd=fracdate(dscd.date_time);
    
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

    [ft]=fracdate(times, 'yyyymmdd_HH:MM:SS');
    times=datetime(times, 'inputformat', 'yyyyMMdd_HH:mm:SS')';
    

    %% kick out bad profiles

    ind=find(prof(1,:)==-9999 | prof(1,:)==-11);
    if ~isempty(ind)
        ft(ind)=[];
        times(ind)=[];
    end

    prof(prof==-9999 | prof==-11)=NaN;
    prof_err(prof_err==-9999 | prof_err==-11)=NaN;

    prof=prof(:,all(~isnan(prof)));
    prof_err=prof_err(:,all(~isnan(prof_err)));

    % same for nd profile
    % columns that failed should be the same as for vmr profile
    prof_nd(prof_nd==-9999 | prof_nd==-11)=NaN;
    prof_nd_err(prof_nd_err==-9999 | prof_nd_err==-11)=NaN;

    prof_nd=prof_nd(:,all(~isnan(prof_nd)));
    prof_nd_err=prof_nd_err(:,all(~isnan(prof_nd_err)));

    %% read retrieval details

    cd('../retrieval_details')

    % get filenames
    tmp = dir('retr_*'); 
    f_info={tmp.name}';

%     if max(size(f_info))~=size(ft,1), error('time mismatch'), end

    % define arrays to hold values
    info=zeros(size(ft,1),3); 
    HHMM=cell(size(ft,1),1);

    info=array2table(info,'VariableNames',{'DOFS','col','col_err'});
    
    for i=1:max(size(f_info))

        % read all data
        tmp=dlmread(f_info{i},':',0,1);

        % save DoF
        info.DOFS(i)=tmp(4);
        
        % save VCD (molec/cm^2) for tg / optical thickness for aer
        info.col(i)=tmp(7);
        
        % save VCD error / optical thickness error
        info.col_err(i)=tmp(8);

        % get time info
        tmp=strsplit(f_info{i},'.');
        % split along underscores
        tmp=strsplit(tmp{1},'_');
        % save time stamp ('HHMM')
        HHMM(i)=tmp(end);
        
    end
    
    HHMM=cell2table(HHMM,'VariableNames',{'HHMM'});
    
    info=[info,HHMM];


    %% save elevations used for each profile
    
    % weighting functions have elevations in the easiest-to-read format
    cd('../weighting_func')
    
    year_tmp=str2double(ymd(1:4));

    % create elevations table: different for pre-2015 data (5/6 deg is lowest)
    if year_tmp>=2015 
        elevs=zeros(size(ft,1),8);
        elevs=array2table(elevs,'VariableNames',...
                          {'el_90','el_30','el_15','el_10','el_5','el_2','el_1','el_m1'});
    else
        if year_tmp==2010
            elevs=zeros(size(ft,1),7);
            elevs=array2table(elevs,'VariableNames',...
                              {'el_90','el_30','el_15','el_10','el_4','el_2','el_1'});
        else
            elevs=zeros(size(ft,1),6);
            elevs=array2table(elevs,'VariableNames',...
                              {'el_90','el_30','el_15','el_10','el_8','el_5_6'});
        end            
    end
    
    % get filenames
    if option=='a'
        tmp = dir('*_ext_20*'); 
    else
        tmp = dir('*wf_*'); 
    end
    
    f_info={tmp.name}';
    
    for i=1:max(size(f_info))

        if option=='a'
            tmp=dlmread(f_info{i},' ',1,1);
            tmp=tmp(1,tmp(1,:)~=0);
        else
            tmp=dlmread(f_info{i},' ',0,1);
            tmp=tmp(1,tmp(1,:)~=0);
        end
    
        % save number of times each elevation angle was used for individual
        % profile
        
        % angles common in all scans
        elevs.el_90(i)=sum(tmp==90+elev_corr);
        elevs.el_30(i)=sum(tmp==30+elev_corr);
        elevs.el_15(i)=sum(tmp==15+elev_corr);
        elevs.el_10(i)=sum(tmp==10+elev_corr);
        
        if year_tmp<2015
           
            if year_tmp==2010
                % different sequence for 2010: 90, 30, 15, 10, 4, 2, 1:
                elevs.el_4(i)=sum(tmp==4+elev_corr);
                elevs.el_2(i)=sum(tmp==2+elev_corr); 
                elevs.el_1(i)=sum(tmp==1+elev_corr);
            else
                % 8 deg always there, some years have 5, some 6 as lowest elev
                % (only one at a time)
                elevs.el_8(i)=sum(tmp==8+elev_corr);
                elevs.el_5_6(i)=sum(tmp==5+elev_corr | tmp==6+elev_corr);
            end
            
        else
            
            elevs.el_5(i)=sum(tmp==5+elev_corr);

            if year_tmp~=2015,
                elevs.el_2(i)=sum(tmp==2+elev_corr); 
            else
                % no 2 deg elevation for 2015, keep table format for consistency
                elevs.el_2(i)=9999; 
            end
            
            elevs.el_1(i)=sum(tmp==1+elev_corr);
            elevs.el_m1(i)=sum(tmp==-1+elev_corr);
            
        end
    end
    
    
    %% read averaging kernels

    cd('../av_kernels')
    
    avk_col=NaN(size(prof));
    avk=NaN(length(alt),length(alt),length(ft));
    
    for i=1:length(ft)
        
        if option=='a'
            tmp=dlmread(['avk_ext_' ymd '_' info.HHMM{i} '.dat'],'',1,4);
        else
            tmp=dlmread(['avk_' ymd '_' info.HHMM{i} '.dat'],'',1,4);
        end
        
        avk_col(:,i)=tmp(:,1);
        avk(:,:,i)=tmp(:,2:end);
        
    end
    
    cd(cur_dir)
    
    %% plot profiles
    
    if show_plot
        figure(1)

        % color plot of profiles during the day
        surf(ft,alt,prof,'EdgeColor','None', 'facecolor', 'interp')
        view(2)
        colormap(jet(300))
        colorbar
        % xlim([ft(1),ft(end)])
        ylim([-0.1,4.1])

        xlabel(['Fractional day, ' ymd(1:4) ' (UTC)'])
        ylabel('Altitude (km)')

        title(ymd)

        % % %% plot measured/modelled DSCDs
        % % 
        % % figure(2)
        % % 
        % % plot(ft_dscd,dscd(:,5),'bo'), hold on
        % % plot(ft_dscd,dscd(:,end),'rx')

        pause(1)
        
    end
    
    %% save results

    f_out=[savedir ymd folder_str '.mat'];

    if ~exist(savedir, 'dir'), mkdir(savedir); end

    save(f_out,'alt','prof','prof_nd','prof_err','prof_nd_err',...
               'ft','times','dscd','ft_dscd','info','elevs','avk','avk_col');

    
end






