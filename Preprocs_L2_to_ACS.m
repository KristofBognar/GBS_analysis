function Preprocs_L2_to_ACS()
%Convert .L2 spectra to .ASC files (extended column format) for QDOAS 3.1
%   Reads all *.L2 files in the current folder and creates new folder with
%   converted .ASC files
%       Parallel preprocs folder if data is already in a preprocs folder
%       New preprocs folder if data is in instr/year directory
%   Original filenames are kept

%% setup

% make list of all L2 (or L3) files
do_L3=true;

if do_L3
    temp = dir('*.L3'); % for UT-GBS 2001, 2003, and 2004; and PEARL-GBS
else
    temp = dir('*.L2'); 
end

f_list = {temp.name}; % cell array of file names

if isempty(f_list), error('Wrong directory, no preprocs files'); end

% split filename and extension 
f_list_1=cell(size(f_list));

for i=1:size(f_list,2)
    temp=strsplit(char(f_list{i}),'.');
    f_list_1{i}=temp{1};
end

% create output folder depending on location of preprocs files
if ~isempty(f_list)
    if strfind(pwd(),'Preprocs')
        if ~exist('../Preprocs_ASC','dir'), mkdir('../Preprocs_ASC'); end
        prefix='../Preprocs_ASC/';
    else 
        if ~exist('Preprocs_ASC','dir'), mkdir('Preprocs_ASC'); end
        prefix='Preprocs_ASC/';
    end
end

% check which instrument based on filename
if f_list{1}(1)=='P'
    instr='PEARL-GBS';
elseif f_list{1}(1)=='U'
    instr='UT-GBS';
end

% check location based on filename
if f_list{1}(2)=='T'
    locname='Toronto Atmospheric Observatory';
    location.longitude = -79.40;
    location.latitude = 43.66;
    location.altitude = 174;
elseif f_list{1}(2)=='E'
    locname='Eureka';    
    location.longitude = -86.416;
    location.latitude = 80.053;
    location.altitude = 610;
elseif f_list{1}(2)=='V'
    locname='Vascony';    
    location.longitude = -106.98;
    location.latitude = 52;
    location.altitude = 516;
elseif f_list{1}(2)=='C'
    locname='Cabauw';    
    location.longitude = 4.898;
    location.latitude = 51.96;
    location.altitude = 200;
elseif f_list{1}(2)=='R'
    locname='Resolute Bay';    
    location.longitude = -94.87;
    location.latitude = 74.68;
    location.altitude = 0;
end

% check if location changes
% if it does, change the file query to pick specific location
loclist=[];
for i=1:length(f_list)
    loclist=[loclist; f_list{1}(2)];
end
if size(unique(loclist),1)>1, error('Multiple locations; select manually'); end

% display checks

% disp(['Reformatting ' instr ' data from ' locname])
% disp(' ')

% measurement modes
mode_vec = {'ZENITH', 'OFFAXIS', 'DIRECT SUN', 'M', 'X', 'HORIZON', 'ALMUCANTAR'};

%% loop over all files
n=0;
for i=1:length(f_list)
    %% read number of records
    num_rec=dlmread(f_list{i}, '\t', [0,0,0,0]);

    %% open output file
    % open file
    if do_L3
        outfile=[prefix f_list_1{i} '_L3.ASC'];
    else
        outfile=[prefix f_list_1{i} '_L2.ASC'];
    end
    
    fid = fopen(outfile, 'w');
    
    %% display progress info
    disp_str=['Reading file ', num2str(i), '/', num2str(length(f_list)), ' (', f_list_1{i}, ')'];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str)+6;    
    
    %% print output file header
    fprintf(fid, '%s\n', ['# Station = ' locname ' (' num2str(location.longitude) ', ' num2str(location.latitude) ')']);
    fprintf(fid, '%s\n', '# Institute = University of Toronto');
    fprintf(fid, '%s\n', '# PI name = Kimberly Strong (strong@atmosp.physics.utoronto.ca)');
    fprintf(fid, '%s\n', ['# Instrument = ' instr]);
    fprintf(fid, '%s\n', '# Size of the detector = 2048');
    fprintf(fid, '%s\n', ['# Total number of records = ' num2str(num_rec)]);
    
    %% loop over all spectra in L2 file
    row_start=2;
    row_end=2056;
    for j=1:num_rec
        % first record starts on line 2
        % each record is 1 + 7 + 2048 = 2056 lines (filled with 0.0 for 2000 pix detector)
        % row of ***, 7 rows of info, and 2048 lines for spectra
        
        % info rows:
        %1 record number
        %2 CCD pixels start (0), CCD pixels end
        %3 mean time, start time, end time (dd mm yyyy HH MM SS, one column each)
        %4 avg solar elev, start solar elev, end solar elev
        %5 shutter, ideal no counts, slit, groove, turret, blaze, centre, 
        %  integration time, no accum, mean TCCD, min TCCD, max TCCD, TBOX 
        %6 meas. mode, viewing elev, viewing az, filter
        %7 location longitude, location latitude
        
        %% read given record
        L2=dlmread(f_list{i}, '\t', [row_start,0,row_end,17]);
        
        % check if reading the right record
        if L2(1,1)~=j, error('Wrong spectra read'), end

        %% assign variabes
        % time
        date=L2(3,1:3); % use mean date
        mean_time=L2(3,4:6);
        start_time=L2(3,10:12);
        end_time=L2(3,16:18);

        % total time (in seconds)
        start_str=[num2str(L2(3,7)) '-' num2str(L2(3,8)) '-' ...
                   num2str(L2(3,9)) ' ' ...
                   num2str(start_time(1)) ':' num2str(start_time(2)) ':' ...
                   num2str(round(start_time(3)))];
        end_str=[num2str(L2(3,13)) '-' num2str(L2(3,14)) '-' ...
                   num2str(L2(3,15)) ' ' ...
                   num2str(end_time(1)) ':' num2str(end_time(2)) ':' ...
                   num2str(round(end_time(3)))];
        total_time=(datenum(end_str,'dd-mm-yyyy HH:MM:SS')-...
                    datenum(start_str,'dd-mm-yyyy HH:MM:SS'))*24*3600;

        % pointing
        v_el=L2(6,2);
        v_az=L2(6,3);
        if v_el<-5 || v_el>90, v_el=9999; end
        if v_az<=0 || v_az>360, v_az=9999; end        
        
        if v_el==90 || v_el==9999 || L2(6,1)==1,
            m_mode=mode_vec{1};
        else
            m_mode=mode_vec{2};
        end
        
        % exposure
        t1int=L2(5,8)/1000; %in seconds
        no_ac=L2(5,9);
        total_exp=t1int*no_ac;

        % sun
        sza=90-L2(4,1);
        saa=get_sol_azim([L2(3,3),L2(3,2),L2(3,1), mean_time], location);
                      
        % spectrum
        L2_spec=L2(8:end,1);
        
        %% new line limits
        row_start=row_start+2056;
        row_end=row_end+2056;
        
        %% write spectrum header
        fprintf(fid,'%s\n', ['Date (DD/MM/YYYY) = ' ...
            sprintf('%02d/%02d/%04d', date(1), date(2), date(3))]);

        fprintf(fid,'%s\n', ['UTC Time (hh:mm:ss) = ' ...
            sprintf('%02d:%02d:%02d', mean_time(1),mean_time(2),round(mean_time(3)))]);

        fprintf(fid,'%s\n', ['UTC Start Time (hh:mm:ss) = ' ...
            sprintf('%02d:%02d:%02d', start_time(1),start_time(2),round(start_time(3)))]);

        fprintf(fid,'%s\n', ['UTC End Time (hh:mm:ss) = ' ...
            sprintf('%02d:%02d:%02d', end_time(1),end_time(2),round(end_time(3)))]);

        fprintf(fid, '%s\n', ['Viewing Elevation Angle (deg) = ' num2str(v_el)]);
        fprintf(fid, '%s\n', ['Viewing Azimuth Angle (deg) = ' num2str(v_az)]);
        fprintf(fid, '%s\n', ['Measurement Type (OFFAXIS/DIRECT SUN/ALMUCANTAR/ZENITH/HORIZON) = ' m_mode]);
        fprintf(fid, '%s\n', ['Exposure time (sec) = ' num2str(t1int)]);
        fprintf(fid, '%s\n', ['Number of scans = ' num2str(no_ac)]);
        fprintf(fid, '%s\n', ['Total Measurement Time (sec) = ' num2str(total_time)]);
        fprintf(fid, '%s\n', ['Total Acquisition Time (sec) = ' num2str(total_exp)]);
        fprintf(fid, '%s\n', ['Solar Zenith Angle (deg) = ' num2str(sza)]);
        fprintf(fid, '%s\n', ['Solar Azimuth Angle (deg) = ' num2str(saa) ...
                                 ' (North=0, East=90)']);

        %% write spectrum
        fprintf(fid, '\t %.6f \n', L2_spec);
        
    end

    % close output file
    fclose(fid);
    
    fprintf('\n');
    fprintf('Done\n');
    
end

function azim = get_sol_azim(date_vec, location)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function sets up the time structure and calls the sun_position
% function to calculate the solar elevation angle
%
% INPUT: date_vec [yyyy, mm, yy, HH, MM, SS]
%        location (structure)
%
% OUTPUT: azim (float azimuth angle of sun)

time.year = date_vec(1);
time.month = date_vec(2);
time.day = date_vec(3);
time.hour =date_vec(4);
time.min = date_vec(5);
time.sec = date_vec(6);
time.UTC = 0;
sun = sun_position(time, location);
azim = sun.azimuth;

end
end

