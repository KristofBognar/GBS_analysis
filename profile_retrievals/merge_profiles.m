function merge_profiles( merge_type, species, year, match_aer, version )
%MERGE_PROFILES merge profile results from individual days/years into a single file
% 
% INPUT:
%   merge_type: 'day': merge daily files for given year
%               'year': merge yearly files for given year range
%               'cindi2': merges daily cindi-2 retrieval results
%                         (automatically cconsiders all versions in
%                         CINDI-2/matlab_files) -- only species input is required 
%               'pandora': merges all data found in folder specified by
%                          'version' variable. Aer and tg data are matched
%                          as well
%
%   species: 'aerosol' for aerosols, 'tracegas' for BrO
%
%   year: individual year for 'day' option; or year range for 'year' option
%
%   match_aer: optional, code will load merged files and remove aerosol
%              profiles that have no matching BrO profile.
%              Prior to 2015, data will be loaded yearly; match_aer==year
%              The other inputs can be left blank for this option
%
%   version: for Eureka: version identifier after eureka_<year> in folder name, empty
%            by default
%            for Pandora: subfolder name (instrument number + input spectra version)
%
% OUTPUT:
%   Saved files with merged profiles
%   Merges _all (all profiles) and _filt (complete scans only) fies separately

% to run for entire period:
%
% for i=2015:2019
%     merge_profiles('day','aerosol',i);
%     merge_profiles('day','tracegas',i);
% end
% merge_profiles('year','aerosol');
% merge_profiles('year','tracegas');
% merge_profiles('','','','1');
%
% @Kristof Bognar, 2018
%
% God this code is terrible  -_-'

if nargin<5
    version=''; 
else
    version=['_' version]; 
end

cur_dir=pwd;
times=[];
info=[];

cd('/home/kristof/work/profile_retrievals/profile_results');

%% filter aerosol files to match BrO profiles
% BrO retrieval fails sometimes (and aer retrieval is missing a few times too)
% Run daily and yearly merges first, this is the last step
% 
% This code OVERWRITES merged files, if they exist
if nargin==4

    disp(['Warning: This option overwrites merged Eureka files']);
    prompt = 'Are you sure you want to contitue? Y/N [N]: ';
    str = input(prompt,'s');
    if ~strcmpi(str,'y')
        return
    end
    
    
    if match_aer==2013 || match_aer==2011 || match_aer==2010  
        % files are not merged from multiple years, just take filtered aer
        % and BrO files for each year, and match the times
        % Creates <tg>_profiles_<year>_all.mat files

        % load BrO times
        load(['eureka_' num2str(match_aer) '/tracegas/profiles_' num2str(match_aer) '_filt.mat'])
        times_bro=times;

        % find common times
        load(['eureka_' num2str(match_aer) '/aerosol/profiles_' num2str(match_aer) '_filt.mat'])

        [~,ind_aer]=setdiff(times,times_bro); % aer times not in BrO
        [~,ind_bro]=setdiff(times_bro,times); % BrO times not in aer

        % filter aerosol data
        filter_data(ind_aer,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,...
                    ['aerosol_profiles_' num2str(match_aer)  '_all.mat'],...
                    alt,dscd,ft_dscd);

        % filter BrO data
        load(['eureka_' num2str(match_aer) '/tracegas/profiles_' num2str(match_aer) '_filt.mat'])
        filter_data(ind_bro,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,...
                    ['tracegas_profiles_' num2str(match_aer) '_all.mat'],...
                    alt,dscd,ft_dscd);
                
    elseif exist('aerosol_profiles_all.mat','file')...
           && exist('tracegas_profiles_all.mat','file')...
           && ~exist('matched_times.txt','file')...
        % Need files that contain merged aerosol and BrO files for
        % 2015-2019
        % These files are loaded and overwritten (times when one
        % retrieval failed are removed)
    
        % load BrO times
        load('tracegas_profiles_all.mat')
        times_bro=times;

        % find common times
        load('aerosol_profiles_all.mat')

        [~,ind_aer]=setdiff(times,times_bro); % aer times not in BrO
        [~,ind_bro]=setdiff(times_bro,times); % BrO times not in aer

        % filter aerosol data
        filter_data(ind_aer,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,'aerosol_profiles_all.mat',...
                    alt,dscd,ft_dscd);

        % filter BrO data
        load('tracegas_profiles_all.mat')
        filter_data(ind_bro,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,'tracegas_profiles_all.mat',...
                    alt,dscd,ft_dscd);

        %%% repeat for _filt files

        % load BrO times
        load('tracegas_profiles_filt_all.mat')
        times_bro=times;

        % find common times
        load('aerosol_profiles_filt_all.mat')

        [~,ind_aer]=setdiff(times,times_bro); % aer times not in BrO
        [~,ind_bro]=setdiff(times_bro,times); % BrO times not in aer

        % filter aerosol data
        filter_data(ind_aer,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,'aerosol_profiles_filt_all.mat',...
                    alt,dscd,ft_dscd);

        % filter BrO data
        load('tracegas_profiles_filt_all.mat')
        filter_data(ind_bro,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,'tracegas_profiles_filt_all.mat',...
                    alt,dscd,ft_dscd);

        % stop here
        fclose(fopen('matched_times.txt', 'w'));
        
    else
        error('Aer and BrO times matched already')
    end

    return

end

%% merge daily or yearly files (Eureka), or CINDI-2 or Pandora profiles
switch merge_type
    
    case 'day' % merge daily files for given year
        
        % change directory
        try
            cd(['eureka_' num2str(year) version '/' species '/' ]);
        catch
            error('Must provide vaild year [merge_profiles(day_or_year,species,years)]');
        end
        
        % get list of files
        tmp=dir('20*.mat');
        f_list={tmp.name}';
        
        % merge files, and create filtered file as well
        merge_and_save(f_list, ['profiles_' num2str(year)], true);
    
    case 'year' % merge yearly files for given year range
        
        % default year range
        if nargin==2
            year=2015:2019; 
            try delete('matched_times.txt'), end
        else
            error('Merge similar years only: use 2015-2019')
        end
        
        f_list_all=cell(length(year),1);
        f_list_filt=cell(length(year),1);
        
        % generate file names
        for i=1:length(year)
            
            f_list_all{i}=['eureka_' num2str(year(i)) '/' species...
                           '/profiles_' num2str(year(i)) '_all.mat'];
            f_list_filt{i}=['eureka_' num2str(year(i)) '/' species...
                           '/profiles_' num2str(year(i)) '_filt.mat'];
            
        end
        
        % merge files
        merge_and_save(f_list_all, [species '_profiles'], false);
        merge_and_save(f_list_filt, [species '_profiles_filt'], false);
        
    case 'cindi2' % merge daily files for given year
        
        % change directory
        cd('CINDI-2/matlab_files/');

        % get list of folders
        tmp=dir('*');
        d_list={tmp.name};
        d_list(1:2)=[];
        
        for d_name=d_list
        
            cd([d_name{1} '/' species '/' ])
            
            % get list of files
            tmp=dir('20*.mat');
            f_list={tmp.name}';

            % merge files, and create filtered file as well
            merge_and_save(f_list, d_name{1}, false);
    
            cd('../../')
            
        end
        
    case 'pandora'
        
        disp('Procesing both aerosol and tracegas data')
        cd(['Pandora/' version(2:end)]);
        pwd_tmp=pwd;

        for species_tmp={'aerosol','tracegas'}
            % change directory
            try
                cd([pwd_tmp '/' species_tmp{1} '/' ]);
            catch
                error('Must provide vaild version [merge_profiles(day_or_year,species,'','',version)]');
            end

            % get list of files
            tmp=dir('20*.mat');
            f_list={tmp.name}';

            % merge files, and create filtered file as well
            merge_and_save(f_list, [species_tmp{1} '_profiles'], true, 'pandora');
        end
        
        %%% match aer and tg profile times
        cd(pwd_tmp)
        % tg
        load('tracegas/tracegas_profiles_all.mat')
        times_tg=times;

        % find common times
        load('aerosol/aerosol_profiles_all.mat')

        [~,ind_aer]=setdiff(times,times_tg); % aer times not in BrO
        [~,ind_tg]=setdiff(times_tg,times); % BrO times not in aer

        % filter aerosol data
        filter_data(ind_aer,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,'aerosol_profiles_all.mat',...
                    alt,dscd,ft_dscd,1);

        % filter tg data
        load('tracegas/tracegas_profiles_all.mat')
        filter_data(ind_tg,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,'tracegas_profiles_all.mat',...
                    alt,dscd,ft_dscd,1);

        %%% repeat for _filt files
        % tg
        load('tracegas/tracegas_profiles_filt.mat')
        times_tg=times;

        % find common times
        load('aerosol/aerosol_profiles_filt.mat')

        [~,ind_aer]=setdiff(times,times_tg); % aer times not in BrO
        [~,ind_tg]=setdiff(times_tg,times); % BrO times not in aer

        % filter aerosol data
        filter_data(ind_aer,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,'aerosol_profiles_filt.mat',...
                    alt,dscd,ft_dscd,1);

        % filter tg data
        load('tracegas/tracegas_profiles_filt.mat')
        filter_data(ind_tg,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,'tracegas_profiles_filt.mat',...
                    alt,dscd,ft_dscd,1);
        
        
end

cd(cur_dir);

end

%% functions that do the work
function merge_and_save(f_list, save_name, filter, filter_type)

    % f_list: full path to files (from current directory)
    % save_name: full path to saved file
    % filter: if true, remove incomplete scans and save separate filtered file
    % default filter option is for eureka PGBS measurements, can also
    % specify pandora instruments
    if nargin==3, filter_type='eureka'; end
    
    % also functions, have to initialize
    info=[];
    times=[];
    
    dscd_all=[];
    times_all=[];
    ft_all=[];
    ft_dscd_all=[];
    info_all=[];
    elevs_all=[];
    avk_all=[];
    avk_col_all=[];

    prof_all=[];
    prof_err_all=[];
    prof_nd_all=[];
    prof_nd_err_all=[];

    for file=f_list'

        try
            load(file{1})
        catch
            error(['File ''' file{1} ''' doesn''t exist']);
        end
            
        if size(info,1)~=size(prof,2), error('Retrieval info missing'); end

        dscd_all=[dscd_all; dscd];
        times_all=[times_all; times];
        ft_all=[ft_all; ft];
        ft_dscd_all=[ft_dscd_all; ft_dscd];
        info_all=[info_all; info];
        elevs_all=[elevs_all; elevs];
        avk_col_all=[avk_col_all,avk_col];

        avk_all=cat(3,avk_all,avk);

        prof_all=[prof_all, prof];
        prof_err_all=[prof_err_all, prof_err];
        prof_nd_all=[prof_nd_all, prof_nd];
        prof_nd_err_all=[prof_nd_err_all, prof_nd_err];


    end

    dscd=dscd_all;
    times=times_all;
    ft=ft_all;
    ft_dscd=ft_dscd_all;
    info=info_all;
    elevs=elevs_all;

    avk_col=avk_col_all;
    avk=avk_all;

    prof=prof_all;
    prof_err=prof_err_all;
    prof_nd=prof_nd_all;
    prof_nd_err=prof_nd_err_all;
    
    save([save_name '_all.mat'],...
         'alt','prof','prof_nd','prof_err','prof_nd_err',...
         'ft','times','dscd','ft_dscd','info','elevs','avk','avk_col');

    if filter

        %%% old version, not extra filter
%         if length(unique(times.Year))==1 % one year only, original filtering works
%             % filter data 
%             if year<=2015 % only include full scans for 2015 and before, except 2010
%                 if year==2010
%                     ind_filt=find(sum(elevs{:,:}==0,2)>1 | elevs.el_90==0);
%                 else
%                     ind_filt=find(sum(elevs{:,:}==0,2)>0);
%                 end
%             else % 2016 and later: allow partial scans with one elev missing(filter by DOFS later)
%                  % must still have 90 deg
%                 ind_filt=find(sum(elevs{:,:}==0,2)>1 | elevs.el_90==0);
%             end
%         
%         else 
        
        % might have multiple years, have to filter separately for each year
        ind_filt=[];
        if strcmp(filter_type,'eureka')
            for yr_tmp=unique(times.Year)'

                if yr_tmp==2010 || yr_tmp>2015
                    % 2016 and later: allow partial scans with one elev
                    % missing (filter by DOFS later)
                    % must still have 90 deg
                    ind_filt=[ind_filt; find((sum(elevs{:,:}==0,2)>1 | elevs.el_90==0) &...
                                             times.Year==yr_tmp)];
                else % only include full scans for 2015 and before, except 2010
                    ind_filt=[ind_filt; find(sum(elevs{:,:}==0,2)>0 & times.Year==yr_tmp)];
                end

            end
            
        elseif strcmp(filter_type,'pandora')
            
            % down-up scans: would contain two of every angle above 1 deg,
            % and one 1deg angle -- exclude scans that miss any elev, or
            % have multiple down-up scans (more than 2 of any angle)
            
            % need to differentiate short and long scans
            % short scans miss most angles, check random angle
            if sum(elevs.el_3)>0 % long scans
                ind_filt=find( sum(elevs{:,:}==0,2)>0 | sum(elevs{:,:}>2,2)>0);
            elseif sum(elevs.el_3)==0 % short scans
                elevs_tmp=elevs(:,{'el_90','el_30','el_15','el_2','el_1'});
                ind_filt=find( sum(elevs_tmp{:,:}==0,2)>0 | sum(elevs_tmp{:,:}>2,2)>0);
            end
        end
        
        filter_data(ind_filt,times,ft,info,elevs,avk_col,avk,...
                    prof,prof_err,prof_nd,prof_nd_err,[save_name '_filt.mat'],...
                    alt,dscd,ft_dscd);
                
    end
     
end

function filter_data(ind_filt,times,ft,info,elevs,avk_col,avk,...
                     prof,prof_err,prof_nd,prof_nd_err,save_name,...
                     alt,dscd,ft_dscd,pan_for_python)

    if nargin==15, pan_for_python=0; end
        
    times(ind_filt)=[];
    ft(ind_filt)=[];
    info(ind_filt,:)=[];
    elevs(ind_filt,:)=[];

    avk_col(:,ind_filt)=[];
    avk(:,:,ind_filt)=[];

    prof(:,ind_filt)=[];
    prof_err(:,ind_filt)=[];
    prof_nd(:,ind_filt)=[];
    prof_nd_err(:,ind_filt)=[];
    
    save(save_name,...
         'alt','prof','prof_nd','prof_err','prof_nd_err',...
         'ft','times','dscd','ft_dscd','info','elevs','avk','avk_col');

     if pan_for_python
        dofs=info.DOFS;
        col=info.col;
        col_err=info.col_err;
        elevs=table2array(elevs);
        dscd=table2array(dscd(:,2:end));
        dscd_head={'SZA','elev','rel_azim','wl','meas','err_meas','retr'};
        
        save([save_name(1:end-4) '_python.mat'],...
             'alt','prof','prof_nd','prof_err','prof_nd_err',...
             'ft','times','dscd','dscd_head','ft_dscd','dofs','col','col_err','elevs','avk','avk_col');
    end
   
end
