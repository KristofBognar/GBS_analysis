function scan_times=insert_90_old( table_in, savedir )
%INSERT_90 Insert 90 deg measurements into QDOAS output files for MAX-DOAS data
%
%
%
% OUTDATED: code is based on splitting dSCDs by local day, then option to
% split by UTC day was tacked on -- not a roboust solution, misses partial
% days and unnecessarily complicated. Use newer version
%
%
%
% Reads MAX-DOAS table and breaks it up by local day (using specNo).
%
% 90 deg placeholders are inseted using interpolated properties from -1 deg and
% next 30 deg measurements
%       - first and last line are extrapolated
%       - partial final scans (ending on 2 deg or lower) are accepted
%       - no correction for missing elev angles
%
% AFTER sorting by day and addition of 90 deg placeholder lines, the scan
% across UTC midnight is deleted and anything on the next day is removed
% from current file (saved for next day)
%       - once sorting and 90 deg lines are done for the next day,
%         past-midnight data is inserted at the top (gap in the morning
%         would mess up time interpolation)
%
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

%% format date/time fields properly
table_in.Timehhmmss.Format='HH:mm:ss';
table_in.DateDDMMYYYY.Format='dd/MM/yyyy';

%% get rid of missing times
ind_nat=find(isnat(table_in.Timehhmmss));
ind_ok=find(~isnat(table_in.Timehhmmss));

table_in.Timehhmmss(ind_nat)=interp1(table_in.Fractionalday(ind_ok),...
                                     table_in.Timehhmmss(ind_ok),...
                                     table_in.Fractionalday(ind_nat));

if sum(isnat(table_in.Timehhmmss)), error('Could not fill all NaTs'); end

scan_times=table;

%% find days

% split into days based on spectrum number (counting restarts from 1 or 2
% for new day)
ind_newday=find(table_in.SpecNo(1:end-1)-table_in.SpecNo(2:end) >0);
% add the first and last day
ind_newday=[1;ind_newday+1; size(table_in,1)];

% get day numbers (Jan 1, 00:00 = 1)
days=floor(table_in.Fractionalday(ind_newday(1:end-1)));

% for breaking up days by UTC as well
past_midnight=[];

%% loop over each day
for i=1:length(days)

    % create daily table
    table_day=table_in(ind_newday(i):ind_newday(i+1)-1,:);
    
    %% first 90
    % check first measurement
    if table_day.Elevviewingangle(1)~=30

        % if not 30 deg, then find actual start of sequence (sometimes
        % there's a -1 deg measurement first)
        ind_start=find(table_day.Elevviewingangle==30);
        ind_start=ind_start(1);
        
%         % throw error if measurements are missing, deal with it if actually happens
%         if ind_start>2, error('Bad first scan'); end

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
    
    table_day(1,14:end)=num2cell(zeros(size(table_day(1,14:end))));
    
    
    %% mid-file 90
    % find places where spectrum number skips 1, that's where 90 deg
    % placeholders should go
    ind90=find(table_day.SpecNo(1:end-1)-table_day.SpecNo(2:end)~=-1);
    ind90=ind90+1;
    
    % check to see whether other elevations were skipped? (unnecessary 90
    % are inserted, doesn't really impact retrievals
    
    % insert 90 deg dummies
    for j=1:length(ind90)
        
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
        table_day(ind90(j),14:end)=num2cell(zeros(size(table_day(ind90(j),14:end))));
        
        % account for added row
        ind90=ind90+1;
    end
    
    
    %% last 90
    % check last measurement (accept partial scans as well)
    if ~any([2,1,0,-1]==table_day.Elevviewingangle(end))
        
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

        table_day(end,14:end)=num2cell(zeros(size(table_day(1,14:end))));
    
    end
    
    %% deal with measurements past UTC midnight
    % append to current file if anything is carried over from previous day
    if ~isempty(past_midnight)
        
        % check for missing day
        if floor(past_midnight.Fractionalday)~=days(i)
            
            % no data for the rest of the day: write what we have to file
            scan_times=write_daily_files(past_midnight,savedir,scan_times);
            
        else
            
            % add to current daily file
            table_day=[past_midnight;table_day];
            
        end
        
        % reset temporary variable
        past_midnight=[];
         
    end
    
    % find last 90 deg line before midnight
    last90=find(floor(table_day.Fractionalday)==days(i) & ...
                table_day.Elevviewingangle==90);
    last90=last90(end);
    
    % if not the last measurement:
    if last90~=size(table_day,1)
    
        % find first 90 deg line after midnight
        first90=find(floor(table_day.Fractionalday)==days(i)+1 & ...
                     table_day.Elevviewingangle==90);
        first90=first90(1);

        % save complete scans past midnight for next day
        if first90~=size(table_day,1)
            
            past_midnight=table_day(first90:end,:);

            % sometimes frac time and date don't match (interp fails when
            % the two scans are on either side of midnight)
            if past_midnight.DateDDMMYYYY(1)~=past_midnight.DateDDMMYYYY(2)
                
                past_midnight.DateDDMMYYYY(1)=past_midnight.DateDDMMYYYY(2);
                past_midnight.Timehhmmss(1)=past_midnight.Timehhmmss(1)-hours(12);
                past_midnight.Fractionaltime(1)=past_midnight.Fractionaltime(1)-12;
                
            end
            
        end
            
        % remove stuff past midnight (including scan across midnight)
        table_day(last90+1:end,:)=[];
        
    end
    
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
    tmp_date=table_day.DateDDMMYYYY+timeofday(table_day.Timehhmmss);
    tmp_date.Format='dd/MM/uuuu HH:mm:ss';
    
    tmp=table;
    tmp.date=tmp_date;
    
    tmp(table_day.Elevviewingangle~=90,:)=[];
    
    tmp.diff=[diff(tmp.date);-100];
    tmp.diff(tmp.diff<0)=NaN;
    
    scan_times=[scan_times;tmp];
end
