function create_retrieval_files()
% Create daily retrieval settings file for aerosol and BrO retrievals
%
% dSCD files created by insert_90.m
%
% Enter start-end times and intervals manually, so files are checked over before the
% retrievals; there are too many exceptions and issues for automating the process
%   
% Aerosol a priori is fixed
% BrO a priori determined from BrO dSCDs and T inversion
%       Surface conc.: 5 when 5deg dscd consistently above 1e14
%       Scale height: from get_sonde_PT.m (check T profile plots as well)
% 
% retrieval performed on Windows side; separate code (python) modifies
% yealy control file with sonde filenames, and copies appropriate retrieval
% file to xx_retrieval.inp
%
% @Kristof Bognar, sometime in 2018
%
%% setup

aer=0; % 1 for aerosol, 0 for BrO
year=2019;

use_uniform_BrO_surf_conc=5; % 0: use original (1 or 5 pptv)
                             % 1: set all surf. conc. to 1 pptv
                             % 5: set all surf. conc. to 5 pptv
                             % change_folder setting takes precedence!

% change number of aerosol iteration steps? 
iter_step=[]; % leave at 5, as in template
% iter_step=10; % increase to given number

% change output folder
change_folder=1;

if change_folder
%     out_folder_name='aer_10iter';
%     out_folder_name='surf_ppt_5to1'; % resets all 5 ppt suf conc to 1
    out_folder_name='surf_ppt_1to5'; % resets all 1 ppt suf conc to 5
else
    out_folder_name='';
    if ~isempty(iter_step), warning('Aer iterations changed, should use different folder'); end
end
    
% start/stop times and BrO a priori are selected manually, and saved
[dates,daily_times,apriori_BrO] = variable_init(year);

%% read template

if aer
    fid=fopen(['/media/kristof/Windows7_OS/SCIATRAN2/AEROSOL_RETRIEVAL_v-1-2/',...
                   'IDL_execute/aerosol_retrieval_template_general.inp'],'r');
else
    fid=fopen(['/media/kristof/Windows7_OS/SCIATRAN2/TRACEGAS_RETRIEVAL_v-1-2/',...
                   'IDL_execute/tracegas_retrieval_template_general.inp'],'r');
end

i = 1;
tline = fgetl(fid);
infile{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    infile{i} = tline;
end
fclose(fid);

%% loop over desired dates
for i=1:length(dates)
    
    
    outfile=infile;
    
    if ~strcmp(daily_times{i,4},datestr(dates(i),'mmdd'))
        error([datestr(dates(i),'mmm dd') ' missing from daily_times cell'])
    end
    
    % special cases:
    % recreate retrieval files only for days when BrO surf conc was set to
    % 5 ppt, ppt, skip other days
    if strcmp(out_folder_name,'surf_ppt_5to1') && ~aer
        
        if i==1, disp('WARNING: all BrO a priori surf. conc. are set to 1 pptv'), end
        
        if str2double(apriori_BrO{i,1})==5
            apriori_BrO{i,1}='1';
        else
            continue
        end
    end
    
    % recreate retrieval files only for days when BrO surf conc was set to
    % 1 ppt, ppt, skip other days
    if strcmp(out_folder_name,'surf_ppt_1to5') && ~aer
        
        if i==1, disp('WARNING: all BrO a priori surf. conc. are set to 5 pptv'), end
        
        if str2double(apriori_BrO{i,1})==1
            apriori_BrO{i,1}='5';
        else
            continue
        end
    end

    % regular files: set BrO a priori surf conc to required value
    if ~change_folder && use_uniform_BrO_surf_conc > 0 && ~aer
        apriori_BrO{i,1}=num2str(use_uniform_BrO_surf_conc);
    end
    
    %% modify relevant lines

    %%% look for field description for dates/times, since placeholder and field
    %%% header repeat
    
    % start date
    tmp=find_in_file(outfile,'Start day of considered time period');
    outfile{tmp+2}=datestr(dates(i),'dd/mm/yyyy');

    % stop date
    tmp=find_in_file(outfile,'Stop day of considered time period');
    outfile{tmp+2}=datestr(dates(i),'dd/mm/yyyy');
    
    % start time
    tmp=find_in_file(outfile,'Start time for the considered day');
    outfile{tmp+2}=daily_times{i,1};

    % stop time
    tmp=find_in_file(outfile,'Stop time for the considered day');
    outfile{tmp+2}=daily_times{i,2};

    % time interval
    tmp=find_in_file(outfile,'Interval length for the temporal retrieval sequences');
    outfile{tmp+2}=daily_times{i,3};
    
    %%% output folder
    if change_folder
        tmp=find_in_file(outfile,'Folder for the output (will be generated if necessary)');
        outfile{tmp+2}=out_folder_name;
    end
    
    %%% look for placeholder for dSCD and sonde files
    
    % dSCD input file
    tmp=find_in_file(outfile,'DSCD_xxxx_xx_xx.dat');
    outfile{tmp}=['DSCD_' datestr(dates(i),'yyyy_mm_dd') '.dat'];

    % radiosonde input file (indices specific to path!!)
    tmp=find_in_file(outfile,'radiosonde_xxxxxx.dat');
    outfile{tmp}(end-9:end-4)=datestr(dates(i),'yymmdd');
    
    % BrO specific stuff
    if ~aer
        
        if ~strcmp(apriori_BrO{i,3},datestr(dates(i),'mmdd'))
            error([datestr(dates(i),'mmm dd') ' missing from apriori_BrO cell'])
        end
        
        % surface concentration in ppm
        tmp=find_in_file(outfile,'First parameter for the a priori trace gas profile');
        outfile{tmp+2}=[apriori_BrO{i,1} 'e-6'];
        
        % profile scale height in km
        tmp=find_in_file(outfile,'Second parameter for the a priori trace gas profile');
        outfile{tmp+2}=apriori_BrO{i,2};
        
        % directory with aerosol results
        if change_folder
            tmp=find_in_file(outfile,'Path to the aerosol retrieval results');
            outfile{tmp+3}=['C:\SCIATRAN2\AEROSOL_RETRIEVAL_v-1-2\Campaign\'...
                            out_folder_name '\'];
        end
        
        
    elseif ~isempty(iter_step)
        
        % change number of iteration steps in aerosol retrieval
        tmp=find_in_file(outfile,'Maximum number of iteration steps');
        outfile{tmp+2}=[num2str(iter_step)];
        
    end
    
    % remove 2deg from list of elevation angles from 2015 data (different
    % scanning sequence for that year)
    if year==2015
        tmp=find_in_file(outfile,'# Elevation angle values to be considered in the retrieval');
        outfile{tmp+2}=['90, 30, 15, 10, 5, 1, -1'];

    % change list of elevation angles and azimuth angle for 2013, 2011, 2010 data
    % (different scanning sequence and tracker location before 2015)
    elseif year==2013
        tmp=find_in_file(outfile,'# Elevation angle values to be considered in the retrieval');
        outfile{tmp+2}=['90, 30, 15, 10, 8, 5'];

        tmp=find_in_file(outfile,'# The single azimuth value when the azimuth is not present');
        outfile{tmp+2}=['35'];
        
        tmp=find_in_file(outfile,'# When azimuth angle selection is performed, select the azimuth angle');
        outfile{tmp+2}=['35'];
    
    elseif year==2011
        tmp=find_in_file(outfile,'# Elevation angle values to be considered in the retrieval');
        outfile{tmp+2}=['90, 30, 15, 10, 8, 6'];

        tmp=find_in_file(outfile,'# The single azimuth value when the azimuth is not present');
        outfile{tmp+2}=['35'];
        
        tmp=find_in_file(outfile,'# When azimuth angle selection is performed, select the azimuth angle');
        outfile{tmp+2}=['35'];
        
        % pointing was off in 2011
        tmp=find_in_file(outfile,'# Value for elevation angle correction of viewing direction');
        outfile{tmp+2}=['1.0'];
        
    elseif year==2010
        tmp=find_in_file(outfile,'# Elevation angle values to be considered in the retrieval');
        outfile{tmp+2}=['90, 30, 15, 10, 4, 2, 1'];

        tmp=find_in_file(outfile,'# The single azimuth value when the azimuth is not present');
        outfile{tmp+2}=['300'];
        
        tmp=find_in_file(outfile,'# When azimuth angle selection is performed, select the azimuth angle');
        outfile{tmp+2}=['300'];
        
        % pointing might be off in 2010
        tmp=find_in_file(outfile,'# Value for elevation angle correction of viewing direction');
        outfile{tmp+2}=['1.0'];
        
    end
    
    %% write file
    
    fpath=['/home/kristof/Drive/data_for_profile_retrieval/input_files/Eureka_'...
           num2str(year) out_folder_name '/'];
       
    if ~exist(fpath,'dir'), mkdir(fpath); end
    
    if aer
        fname=['aerosol_retrieval_' datestr(dates(i),'yymmdd') '.inp'];
    else
        fname=['tracegas_retrieval_' datestr(dates(i),'yymmdd') '.inp'];
    end
    
    file_tmp = fopen([fpath fname],'w');
    fprintf(file_tmp,'%s\n',outfile{:});
    fclose(file_tmp);        
        

end

end



function [dates,daily_times,apriori_BrO] = variable_init(year)

    %%% start and end times + time interval (match manually to times in dSCD file)
    % (start,stop,interval, mmdd, comments), and
    %%% BrO a priori
    % determined by visual insection of dSCDs and T inversion (Zhao et al., 2016)
    % surf VMR (x1e-6), scake height (km), mmdd
        
    if year==2019
        
        % day range 
        doy_range=[64:151]; 

        % remove missing days
        missing=[112,118];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);
        
        % times and apriori info        
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2019.mat')
       
    end
    
    if year==2018
        %% dates

        % day range 
        doy_range=[64:151]; 

        % remove missing days
        doy_range(doy_range==133)=[]; 

        % dates
        dates=ft_to_date(doy_range-1,year);

        % times and apriori info        
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2018.mat')
                 
    elseif year==2017
        %% dates

        % day range 
        doy_range=[66:94]; 

        % remove missing days
%         doy_range(doy_range==133)=[]; 

        % dates
        dates=ft_to_date(doy_range-1,year);
        
        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2017.mat')       
   
    elseif year==2016
        %% dates

        % day range 
        doy_range=[66:132]; 

        % remove missing days
        missing=[73:75,84:91];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);

        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2016.mat')       

    elseif year==2015
        %% dates

        % day range 
        doy_range=[62:151]; 

        % remove missing days
        missing=[63,87,98,99,106,123,145];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);

        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2015.mat')       
      
    elseif year==2013
        %% dates

        % day range 
        doy_range=[69:113]; 

        % remove missing days
        doy_range(doy_range==93)=[]; 

        % dates
        dates=ft_to_date(doy_range-1,year);
                
        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2013.mat')       

    elseif year==2011
        %% dates

        % day range 
        doy_range=[70:100]; 

        % remove missing days
        missing=[83,84,96,98];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);
                
        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2011.mat')       
        
    elseif year==2010
        %% dates

        % day range 
        %%% there are more measurements; only process part of the dataset
        doy_range=[100:140]; 

        % remove days that are not needed
        missing=[107:122,130];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);
                
        % times and apriori info
        % used 1h profile time window for most days, since scans are very
        % long for some reason (there are lots of gaps too)
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2010.mat')       
        

    end
end






