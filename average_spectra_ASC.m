% Script to average PEARL-GBS spectra taken from L2.ASC files
% Number of spectra to average is hard coded
% Spectra with SZA>86 are not averaged, since MAX-DOAS mode means ZS
%   spectra are spaced too far apart
%
%   INPUT: all L2.ASC files in current folder
%   OUTPUT: averaged spectra in corresponding L3.ASC files
%
% created by Kristof Bognar, September 2017


%% number of spectra to average
num_avg=5;

%% create list of L2 files
% only average zenith measurements
tmp = dir('PE0Z*L2.ASC'); 
f_list = {tmp.name}; % cell array of file names

if isempty(f_list), error('Must be in Preprocs_ASC directory'); end

%% print out instrument and year the setup is optimized for

disp('Script is set up for PEARL-GBS UV data from 2010');
disp('-----');
disp('SEVERAL DAYS ARE AUTOMATICALLY REMOVED!');
disp('-----');
prompt = 'Press Return to continue: ';
str = input(prompt,'s');
if ~isempty(str), return, end

% optional: remove days that were excuded from originl analysis
% list of days to exclude, based on analysis of L2 files
exclude_2010=[79,90,107,108,120,121,122,141,142,143,144,145,146,147,148,162,170,...
              176,179,180,242,243,265,297];

% remove bad days from name list
for i=1:length(exclude_2010)
    
    str=num2str(exclude_2010(i));
    IndexC = strfind(f_list, str);
    ind_day = find(not(cellfun('isempty', IndexC)));
    
    if ~isempty(ind_day) 
        ind_day=ind_day(1);
        f_list(ind_day)=[]; 
    end
    
end

%% create L3 file names and output directory
% split filename and extension 
f_list_out=cell(size(f_list));

for i=1:size(f_list,2)
    tmp=strsplit(char(f_list{i}),'_');
    f_list_out{i}=[tmp{1} '_' tmp{2} '_' tmp{3} '_L3.ASC'];
end

% create output directory
out_dir=['L3_' num2str(num_avg) 'spec'];
if ~exist(out_dir,'dir'), mkdir(out_dir); end



%% loop over all files
n=0;
for file=1:length(f_list)

    % display progress info
    disp_str=['Reading file ',num2str(file),'/',num2str(length(f_list)),' (',f_list{file},')'];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str);    
    
    
    % open files
    fid_in=fopen(f_list{file},'r');

    fid_out=fopen([out_dir '/'  f_list_out{file}],'w');

    % read file info
    comments=cell(6,1);
    for i=1:6
        comments{i}=fgetl(fid_in);
        fprintf(fid_out, '%s\n', comments{i});
    end

    % number of records
    tmp=strsplit(comments{end},'=');
    num_rec=str2num(tmp{2});

    % variables to hold spectra and details
    header=cell(13,1);
    spec=NaN(2048,1);
    % sza=zeros(num_rec);

    %% loop over spectra in the file
    % each record is 13 header lines + 2048 spectra lines
    reinit=true;
    end_avg=false;
    num_rec_new=0;
    for i=1:num_rec

        % read header
        for j=1:13
            header{j}=fgetl(fid_in);
        end

        % read spectra and check if it worked
        % (textscan advances line count for fgetl as well)
        specline=6+(i-1)*2061+13+1;
        spec=cell2mat(textscan(fid_in,'%.6f\n',2048)); 
        
        if isempty(spec) || length(spec)~=2048
            error(['Spectrum not read, line ' num2str(specline)]);
        end

        % set up variables for averaging -- reinitialize after each averaging
        if reinit
            count=1;

            exp_time=NaN(1,num_avg);
            no_scans=NaN(1,num_avg);
            tot_ac_time=NaN(1,num_avg);
            sza=NaN(1,num_avg);
            saa=NaN(1,num_avg);
            date=cell(1,num_avg);

            spec_arr=NaN(2048,num_avg);
        end

        % get SZA
        tmp=strsplit(header{12},'=');
        sza(count)=str2num(tmp{2});


        if sza(count)<=86
            %% if sza<86 (MAX-DOAS mode), write spectrum into L3 file straight away
            % write header lines without modification
            for j=1:13
                fprintf(fid_out, '%s\n', header{j});
            end
            % write spectrum without modification
            fprintf(fid_out, '\t %.6f \n', spec);
            
            % count number of records in file
            num_rec_new=num_rec_new+1;

            % reinitialize variables for averaging
            reinit=true;
            % if sza has decreased below 86 but there are still some spectra
            % (< num_avg) waiting to be averaged
            if count>1, end_avg=true; end

        elseif sza(count)>86
            %% if sza>86, save properties for averaging

            % save exposure details
            tmp=strsplit(header{8},'=');
            exp_time(count)=str2num(tmp{2});

            tmp=strsplit(header{9},'=');
            no_scans(count)=str2num(tmp{2});

            tmp=strsplit(header{11},'=');
            tot_ac_time(count)=str2num(tmp{2});

            % extract SAA info
            tmp=strsplit(header{13},'=');
            tmp2=strsplit(tmp{2},'(');
            saa(count)=str2num(tmp2{1});

            % save date
            tmp=strsplit(header{1},'=');
            date{count}=tmp{2};

            % save start time if it's the first measurement
            if count==1
                tmp=strsplit(header{3},'=');
                start_time=tmp{2};
            end

            % save end time, and overwrite prev. value (only want to keep the last value)
            tmp=strsplit(header{4},'=');
            end_time=tmp{2};

            % save spectrum
            spec_arr(:,count)=spec;

            reinit=false;
            count=count+1;
        end

        %% if number of spectra reaches num_avg, average spectra and write to file
        if count==num_avg+1

            % check if arrays are full (enough spectra to average) -- if
            % average is cut off, then only use the values present
            if end_avg
                lim=count-1;
            else
                lim=num_avg;
            end

            % calculate new measurement time
            tmp1=datenum([date{1} ' ' start_time],'dd/mm/yyyy HH:MM:SS');
            tmp2=datenum([date{lim} ' ' end_time],'dd/mm/yyyy HH:MM:SS');          

            mean_date=datestr((tmp1+tmp2)/2,'dd/mm/yyyy'); % in case of midnight measurements
            mean_time=datestr((tmp1+tmp2)/2,'HH:MM:SS');

            % total measurement time in seconds
            meas_time=(tmp2-tmp1)*24*3600;

            % new measurement parameters
            mean_exp_time=mean(exp_time(1:lim));
            total_no_scans=sum(no_scans(1:lim));
            sum_tot_ac_time=sum(tot_ac_time(1:lim));
            mean_sza=mean(sza(1:lim));
            mean_saa=circ_mean(saa(1:lim)'*(pi/180))*(180/pi); % in case of mignight measurenemts
            if mean_saa<0, mean_saa=360+mean_saa; end

            % check if data was read/saved properly
            if isnan(mean_exp_time) || isnan(total_no_scans) || ...
               isnan(sum_tot_ac_time) || isnan(mean_sza) || isnan(mean_saa)
                error('Spectra info not recorded properly')
            end

            % print header info
            fprintf(fid_out,'%s\n', ['Date (DD/MM/YYYY) = ' mean_date]);
            fprintf(fid_out,'%s\n', ['UTC Time (hh:mm:ss) = ' mean_time]);
            fprintf(fid_out,'%s\n', ['UTC Start Time (hh:mm:ss) =' start_time]);
            fprintf(fid_out,'%s\n', ['UTC End Time (hh:mm:ss) =' end_time]);
            fprintf(fid_out, '%s\n', header{5}); 
            fprintf(fid_out, '%s\n', header{6}); % only Z files, pointing is constant
            fprintf(fid_out, '%s\n', header{7});
            fprintf(fid_out, '%s\n', ['Exposure time (sec) = ' num2str(mean_exp_time)]);
            fprintf(fid_out, '%s\n', ['Number of scans = ' num2str(total_no_scans)]);        
            fprintf(fid_out, '%s\n', ['Total Measurement Time (sec) = ' num2str(meas_time)]);
            fprintf(fid_out, '%s\n', ['Total Acquisition Time (sec) = ' num2str(sum_tot_ac_time)]);
            fprintf(fid_out, '%s\n', ['Solar Zenith Angle (deg) = ' num2str(mean_sza)]);
            fprintf(fid_out, '%s\n', ['Solar Azimuth Angle (deg) = ' num2str(mean_saa) ...
                                     ' (North=0, East=90)']);

            % print averaged spectrum
            fprintf(fid_out, '\t %.6f \n', mean(spec_arr(:,1:lim),2));
            
            % count number of records in file
            num_rec_new=num_rec_new+1;
            
            % reset variables for averaging
            reinit=true;

        end

    end


    fclose(fid_in);
    fclose(fid_out);

    %% replace record number in L3 file
    
    % read in entire file
    fid=fopen([out_dir '/'  f_list_out{file}],'r');
    f=fread(fid,'*char')';
    fclose(fid);    
    
    % replace row containing record number
    f = strrep(f,['# Total number of records = ' num2str(num_rec)],...
                 ['# Total number of records = ' num2str(num_rec_new)]);
    
    % rewrite file
    fid=fopen([out_dir '/'  f_list_out{file}],'w');
    fprintf(fid,'%s',f);
    fclose(fid);    
    
end

fprintf('\n');
fprintf('Done\n');

