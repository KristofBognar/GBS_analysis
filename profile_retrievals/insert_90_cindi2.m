function insert_90_cindi2( table_in, savedir, uvvis)
%INSERT_90 Insert 90 deg measurements into QDOAS output files for MAX-DOAS data
%   
% INPUT:    table_in: MAX-DOAS dSCD file (straight from QDOAS), imported as a
%               matlab table
%           savedir: directory where output files will be saved
%           uvvis: string for file names, 'uv' or 'vis'
%
% OUTPUT:   Daily dSCD files with dummy 90 deg lines inserted
%
%
% Kristof Bognar, March 2018

%% find days

% split into days based on day number 

% get day numbers (Jan 1, 00:00 = 0)
days=floor(table_in.DOY);

ind_newday=find(days(1:end-1)-days(2:end) <0);
% add the first and last day
ind_newday=[1;ind_newday+1; size(table_in,1)+1];


%% add date and time info

dates=datetime(table_in.DOY+yeartime(2016),'convertfrom','datenum');

% table_in.DateDDMMYYYY=datestr(dates,'dd/mm/yyyy');
% table_in.Timehhmmss=datestr(dates,'HH:MM:SS');

table_in.DateDDMMYYYY=dates;
table_in.DateDDMMYYYY.Format='dd/MM/yyyy';
table_in.Timehhmmss=dates;
table_in.Timehhmmss.Format='HH:mm:ss';


table_in=[table_in(:,1:6), table_in(:,11:12), table_in(:,7:10)];


%% loop over each day
for i=1:length(unique(days))

    % create daily table
    table_day=table_in(ind_newday(i):ind_newday(i+1)-1,:);
    
    %% first 90
    % replicate first line, and modify time so it passes as 90 deg dummy
    table_day=[table_day(1,:); table_day];
    
    table_day.DOY(1)=table_day.DOY(1)-diff(table_day.DOY(2:3));
    table_day.UTC(1)=table_day.UTC(1)-diff(table_day.UTC(2:3));
    table_day.SZA(1)=table_day.SZA(1)-diff(table_day.SZA(2:3));
    table_day.SolarAzimuthAngle(1)=table_day.SolarAzimuthAngle(1)-diff(table_day.SolarAzimuthAngle(2:3));
    
    table_day.Elevviewingangle(1)=90;
    
    table_day.Timehhmmss(1)=table_day.Timehhmmss(1)-diff(table_day.Timehhmmss(2:3));
    
    table_day(1,9:end)=num2cell(zeros(size(table_day(1,9:end))));
    
    
    %% mid-file 90
    % find places where spectrum number skips 1, that's where 90 deg
    % placeholders should go
    ind90=find(table_day.Elevviewingangle(1:end-1)-table_day.Elevviewingangle(2:end)==29);
    ind90=ind90+1;
    
    % insert 90 deg dummies
    for j=1:length(ind90)
        
        table_day=[table_day(1:ind90(j)-1,:);...
                   table_day(ind90(j)-1,:);...
                   table_day(ind90(j):end,:)];
        
% %         % average spectra details (frac time, sza, saa, etc)
% %         table_day(ind90(j),1:4)=num2cell(mean(table_day{ind90(j)-1:2:ind90(j)+1,1:4}));
% %         % average datetime fields
% %         table_day(ind90(j),7:8)=num2cell(mean(table_day{ind90(j)-1:2:ind90(j)+1,7:8}));
% %         % add 90 for elev viewing angle
% %         table_day.Elevviewingangle(ind90(j))=90;
% %         % set dSCD data to 0
% %         table_day(ind90(j),9:end)=num2cell(zeros(size(table_day(ind90(j),9:end))));
        
        table_day.DOY(ind90(j))=table_day.DOY(ind90(j))+diff(table_day.DOY(ind90(j)-2:ind90(j)-1));
        table_day.UTC(ind90(j))=table_day.UTC(ind90(j))+diff(table_day.UTC(ind90(j)-2:ind90(j)-1));
        table_day.SZA(ind90(j))=table_day.SZA(ind90(j))+diff(table_day.SZA(ind90(j)-2:ind90(j)-1));
        table_day.SolarAzimuthAngle(ind90(j))=table_day.SolarAzimuthAngle(ind90(j))+diff(table_day.SolarAzimuthAngle(ind90(j)-2:ind90(j)-1));

        table_day.Elevviewingangle(ind90(j))=90;

        table_day.Timehhmmss(ind90(j))=table_day.Timehhmmss(ind90(j))+diff(table_day.Timehhmmss(ind90(j)-2:ind90(j)-1));

        table_day(ind90(j),9:end)=num2cell(zeros(size(table_day(1,9:end))));



        % account for added row
        ind90=ind90+1;
    end
    
    
    %% last 90
    % check last measurement (accept partial scans as well)
    % replicate last line, and modify time so it passes as 90 deg dummy
    table_day=[table_day; table_day(end,:)];

    table_day.DOY(end)=table_day.DOY(end)+diff(table_day.DOY(end-2:end-1));
    table_day.UTC(end)=table_day.UTC(end)+diff(table_day.UTC(end-2:end-1));
    table_day.SZA(end)=table_day.SZA(end)+diff(table_day.SZA(end-2:end-1));
    table_day.SolarAzimuthAngle(end)=table_day.SolarAzimuthAngle(end)+diff(table_day.SolarAzimuthAngle(end-2:end-1));

    table_day.Elevviewingangle(end)=90;

    table_day.Timehhmmss(end)=table_day.Timehhmmss(end)+diff(table_day.Timehhmmss(end-2:end-1));

    table_day(end,9:end)=num2cell(zeros(size(table_day(1,9:end))));
        
    %% save daily file
    
    yyyy=num2str(year(table_day.DateDDMMYYYY(1)));
    mm=num2str(month(table_day.DateDDMMYYYY(1)));
    dd=num2str(day(table_day.DateDDMMYYYY(1)));
    
    if length(mm)==1, mm=['0' mm]; end
    if length(dd)==1, dd=['0' dd]; end
    
    fname=[savedir 'DSCD_' yyyy '_' mm '_' dd '_' uvvis '_v2.csv'];
    
    writetable(table_day,fname,'Delimiter',',');
    
end

end




