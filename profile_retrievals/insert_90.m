function scan_times=insert_90( table_in, savedir )
%INSERT_90 Insert 90 deg measurements into QDOAS output files for MAX-DOAS data
%
% Reads MAX-DOAS table and breaks it up by UTC day.
%
% 90 deg placeholders are inseted using interpolated properties from lowest elev and
% next 30 deg measurements
%       - first and last line are extrapolated
%       - partial final scans (ending on 2 deg or lower) are accepted
%       - no correction for missing elev angles
%
% Scans across midnight are most likely deleted since retrieval requires
% data for a single day (we have daily sonde and BrO a priori)
%
% Code assumes that scans are in descending order by elevation angle!
%   
% INPUT:    table_in: MAX-DOAS dSCD file (straight from QDOAS), imported as a
%               matlab table
%           savedir: directory where output files will be saved
%
% OUTPUT:   Daily dSCD files with dummy 90 deg lines inserted
%           Scan lengths (90 to next 90)
%
%
% Kristof Bognar, March 2018

%% setup
% find first retrieval column (measurement info columns come first; assume
% that first retrieval window is NO2)
% general info columns might change as people add columns with extra info
retr_start_col=find(strcmp(table_in.Properties.VariableNames,'NO2RMS'));

% make sure format of first few columns is the same (if not, change section
% that adds dummy  90deg lines in middle of file)
%%% Expected columns:
% Spec No
% Year
% Fractional day
% Fractional time
% Scans
% Tint
% SZA
% Solar Azimuth Angle
% Elev. viewing angle
% Azim. viewing angle
% Date (DD/MM/YYYY)
% Time (hh:mm:ss)
% Total Experiment Time (sec)

% lowest few elevations (to check for completed last scans)
low_elevs=unique(table_in.Elevviewingangle)';
low_elevs=low_elevs(1:4);

% gap tolerance
% max accepted time difference between consecutive measurements -- if time
% is greater, it's considered a gap
if table_in.Year(1)==2010
    gap_tolerance=0.35; % 21 minutes, measurements are strangely long in 2010
else
    gap_tolerance=0.2; % 12 minutes
end

%% format date/time fields properly
table_in.Timehhmmss.Format='HH:mm:ss';
table_in.DateDDMMYYYY.Format='dd/MM/yyyy';

table_in.Timehhmmss.Year=table_in.DateDDMMYYYY.Year;
table_in.Timehhmmss.Month=table_in.DateDDMMYYYY.Month;
table_in.Timehhmmss.Day=table_in.DateDDMMYYYY.Day;

%% get rid of missing times
ind_nat=isnat(table_in.Timehhmmss);
ind_ok=~isnat(table_in.Timehhmmss);

table_in.Timehhmmss(ind_nat)=...
         ft_to_date(table_in.Fractionalday(ind_nat)-1,table_in.Year(ind_nat));

if sum(isnat(table_in.Timehhmmss)), error('Could not fill all NaTs'); end

scan_times=table;

%% find days

days=unique(table_in.DateDDMMYYYY);

%% loop over each day
for i=1:length(days)

    % create daily table
    table_day=table_in(table_in.DateDDMMYYYY==days(i),:);
    
    % skip if too few measurements
    if size(table_day,1)<8, continue, end 
    
    %% first 90
    % check first measurement
    if table_day.Elevviewingangle(1)~=30

        % if not 30 deg, then find actual start of sequence (sometimes
        % there's a -1 deg measurement first, or the end of the scan from
        % the pervious day)
        ind_start=find(table_day.Elevviewingangle==30);
        ind_start=ind_start(1);
        
        % check if there are actually measurements after the 30 deg line
        % (sometimes there are not)
        search=1;
        while search
            if diff(table_day.Fractionalday(ind_start:ind_start+1))<gap_tolerance
                % scan actually strarts, exit loop
                search=0; 
            else
                % find start of next scan -- accept any angle
                table_day(ind_start,:)=[];
                ind_start=ind_start+1;
            end
        end

        % discard extra measurement(s)
        table_day(1:ind_start-1,:)=[];

    end
        
    % replicate first line, and modify time so it passes as 90 deg dummy
    table_day=[table_day(1,:); table_day];
    
    table_day.SpecNo(1)=table_day.SpecNo(1)-1;
    table_day.Fractionalday(1)=table_day.Fractionalday(1)-diff(table_day.Fractionalday(2:3));
    table_day.Fractionaltime(1)=table_day.Fractionaltime(1)-diff(table_day.Fractionaltime(2:3));
    table_day.SZA(1)=table_day.SZA(1)-diff(table_day.SZA(2:3));
    table_day.SolarAzimuthAngle(1)=table_day.SolarAzimuthAngle(1)-diff(table_day.SolarAzimuthAngle(2:3));
    
    table_day.Elevviewingangle(1)=90;
    
    table_day.Timehhmmss(1)=table_day.Timehhmmss(1)-diff(table_day.Timehhmmss(2:3));
    
    table_day(1,retr_start_col:end)=num2cell(zeros(size(table_day(1,retr_start_col:end))));
    
    % check if time slipped back to prev. day
    if hour(table_day.Timehhmmss(1)) == 23
        
        table_day.Timehhmmss(1).Hour=0;
        table_day.Timehhmmss(1).Minute=0;
        table_day.Timehhmmss(1).Second=2;
        table_day.Fractionaltime(1)=0.0005556;
        
    end
    
    %% mid-file 90
    % find places where spectrum number skips 1 (or more), that's where 90 deg
    % placeholders should go
    ind90=find(table_day.SpecNo(1:end-1)-table_day.SpecNo(2:end)~=-1);
    ind90=ind90+1;
    
    % add places where 90deg spectrum wasn't taken (no gap in spectrum number -- rare)
    ind30=find(table_day.Elevviewingangle==30);
    ind30(1)=[]; % exclude first 30, already has prior 90deg dummy added
    
    ind90=union(ind90,ind30);
    
    % handle data gaps (measurements past midnight and then in afternoon, or missing scans)
    % take gaps greater than 12 min to indicate missing data or measurement break
    ind90_gap=find(table_day.Fractionaltime(1:end-1)-table_day.Fractionaltime(2:end) < -gap_tolerance);
    ind90_gap=ind90_gap+1;
    
    % make sure last gap (if any) has at least two measurements after it
    if ~isempty(ind90_gap) && (size(table_day,1)-ind90_gap(end)<2) 
        % remove corresponding 90 index, so code doesn't attempt to add 90
        % deg measurement
        ind90(ind90==ind90_gap(end))=[];
        table_day(ind90_gap(end):end,:)=[];
        ind90_gap(end)=[];
        
        search=1;
        while search
            if ~isempty(ind90_gap) && (size(table_day,1)-ind90_gap(end)<2) 
                ind90(ind90==ind90_gap(end))=[];
                table_day(ind90_gap(end):end,:)=[];
                ind90_gap(end)=[];
            else
                search=0;
            end
        end    
    end

    % make sure there's enough measurements at the end of the day to
    % interpolate dummy times
    if length(ind90)>1 && ind90(end)==size(table_day,1) 
        % ind90 indicates start of scan
        % only one measurement in last scan: delete measurement and index
        ind90(end)=[];
        table_day(end,:)=[]; 
    end
    
    
    % insert 90 deg dummies
    for j=1:length(ind90)
        
        % check if there's a gap
        if ismember(ind90(j),ind90_gap)
        
            % if yes, insert two dummies, on each end of the gap
            table_day=[table_day(1:ind90(j)-1,:);...
                       table_day(ind90(j)-1,:);...
                       table_day(ind90(j),:);...
                       table_day(ind90(j):end,:)];

            gap_start=ind90(j);
            gap_end=ind90(j)+1;

            % add dummy at start of the gap
            table_day.SpecNo(gap_start)=table_day.SpecNo(gap_start)+1;
            table_day.Fractionalday(gap_start)=table_day.Fractionalday(gap_start)+diff(table_day.Fractionalday(gap_start-2:gap_start-1));
            table_day.Fractionaltime(gap_start)=table_day.Fractionaltime(gap_start)+diff(table_day.Fractionaltime(gap_start-2:gap_start-1));
            table_day.SZA(gap_start)=table_day.SZA(gap_start)+diff(table_day.SZA(gap_start-2:gap_start-1));
            table_day.SolarAzimuthAngle(gap_start)=table_day.SolarAzimuthAngle(gap_start)+diff(table_day.SolarAzimuthAngle(gap_start-2:gap_start-1));

            table_day.Elevviewingangle(gap_start)=90;

            table_day.Timehhmmss(gap_start)=table_day.Timehhmmss(gap_start)+diff(table_day.Timehhmmss(gap_start-2:gap_start-1));

            table_day(gap_start,retr_start_col:end)=num2cell(zeros(size(table_day(1,retr_start_col:end))));

            % add dummy at the end of the gap
            table_day.SpecNo(gap_end)=table_day.SpecNo(gap_end)-1;
            table_day.Fractionalday(gap_end)=table_day.Fractionalday(gap_end)-diff(table_day.Fractionalday(gap_end+1:gap_end+2));
            table_day.Fractionaltime(gap_end)=table_day.Fractionaltime(gap_end)-diff(table_day.Fractionaltime(gap_end+1:gap_end+2));
            table_day.SZA(gap_end)=table_day.SZA(gap_end)-diff(table_day.SZA(gap_end+1:gap_end+2));
            table_day.SolarAzimuthAngle(gap_end)=table_day.SolarAzimuthAngle(gap_end)-diff(table_day.SolarAzimuthAngle(gap_end+1:gap_end+2));

            table_day.Elevviewingangle(gap_end)=90;

            table_day.Timehhmmss(gap_end)=table_day.Timehhmmss(gap_end)-diff(table_day.Timehhmmss(gap_end+1:gap_end+2));

            table_day(gap_end,retr_start_col:end)=num2cell(zeros(size(table_day(1,retr_start_col:end))));

            % account for added rows
            ind90=ind90+2;
            ind90_gap=ind90_gap+2;
        
        else % if no, interpolate between two scans
            
            table_day=[table_day(1:ind90(j)-1,:);...
                       table_day(ind90(j),:);...
                       table_day(ind90(j):end,:)];

            % average spectra details (frac time, sza, saa, etc)
            table_day(ind90(j),1:8)=num2cell(mean(table_day{ind90(j)-1:2:ind90(j)+1,1:8}));
            % average datetime fields
            table_day(ind90(j),11:12)=num2cell(mean(table_day{ind90(j)-1:2:ind90(j)+1,11:12}));
            % add 90 for elev viewing angle
            table_day.Elevviewingangle(ind90(j))=90;
            % set dSCD data to 0
            table_day(ind90(j),retr_start_col:end)=num2cell(zeros(size(table_day(ind90(j),retr_start_col:end))));

            % account for added row
            ind90=ind90+1;
            ind90_gap=ind90_gap+1;
        
        end
    end
    
    
    %% last 90
    % check last measurement (accept partial scans as well)
    if ~any(low_elevs==table_day.Elevviewingangle(end))
        
        % delete incomplete scan at the end, 90deg line has been created by
        % code above
        table_day(ind90(end):end,:)=[];

    else
        
        % replicate last line, and modify time so it passes as 90 deg dummy
        table_day=[table_day; table_day(end,:)];

        table_day.SpecNo(end)=table_day.SpecNo(end)+1;
        table_day.Fractionalday(end)=table_day.Fractionalday(end)+diff(table_day.Fractionalday(end-2:end-1));
        table_day.Fractionaltime(end)=table_day.Fractionaltime(end)+diff(table_day.Fractionaltime(end-2:end-1));
        table_day.SZA(end)=table_day.SZA(end)+diff(table_day.SZA(end-2:end-1));
        table_day.SolarAzimuthAngle(end)=table_day.SolarAzimuthAngle(end)+diff(table_day.SolarAzimuthAngle(end-2:end-1));

        table_day.Elevviewingangle(end)=90;

        table_day.Timehhmmss(end)=table_day.Timehhmmss(end)+diff(table_day.Timehhmmss(end-2:end-1));

        table_day(end,retr_start_col:end)=num2cell(zeros(size(table_day(1,retr_start_col:end))));
    
    end
    
    % check if time slipped to next day
    if hour(table_day.Timehhmmss(end)) == 0
        
        table_day.Timehhmmss(end).Hour=23;
        table_day.Timehhmmss(end).Minute=59;
        table_day.Timehhmmss(end).Second=58;
        table_day.Fractionaltime(end)=23.9994444;
        
    end
    
    %% Recalculate datetime column
    
    table_day.DateTime=table_day.DateDDMMYYYY+timeofday(table_day.Timehhmmss);
    table_day.DateTime.Format='dd/MM/uuuu HH:mm:ss';
    
    %% Save dSCD file
    
    if ~exist(savedir,'dir'), mkdir(savedir); end
    scan_times=write_daily_files(table_day,savedir,scan_times);
    
end

end


function scan_times=write_daily_files(table_day,savedir,scan_times)
    %% save daily file
    
    yyyy=num2str(year(table_day.DateDDMMYYYY(1)));
    mm=num2str(month(table_day.DateDDMMYYYY(1)));
    dd=num2str(day(table_day.DateDDMMYYYY(1)));
    
    if length(mm)==1, mm=['0' mm]; end
    if length(dd)==1, dd=['0' dd]; end
    
    fname=[savedir 'DSCD_' yyyy '_' mm '_' dd '.dat'];
    
    writetable(table_day,fname,'Delimiter',',');
    
    %% calculate and save length of each scan (90 to next 90)
    tmp=table;
    tmp.date=table_day.DateTime;
    
    tmp(table_day.Elevviewingangle~=90,:)=[];
    
    tmp.diff=[diff(tmp.date);-100];
    tmp.diff(tmp.diff<0)=NaN;
    
    scan_times=[scan_times;tmp];
end
