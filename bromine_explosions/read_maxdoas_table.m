function maxdoas = read_maxdoas_table(path,filename)
%read_maxdoas_table Import maxdoas data from QDOAS .ASC file as a table.
%
% Auto-generated by MATLAB on 2019/02/13 16:31:21
% Modified by Kristof Bognar:
%   Add datetime column
%   Format time column to show time instead of date
%   Save output in same folder as the ASC file, under the same name


%% Initialize variables.
delimiter = '\t';
startRow = 3;
endRow = inf;

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
try
    fileID = fopen([path filename],'r');
catch
    error(['Invalid filename ' path filename]);
end

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

dateFormats = {'dd/MM/yyyy', 'HH:mm:ss'};
dateFormatIndex = 1;
blankDates = cell(1,size(raw,2));
anyBlankDates = false(size(raw,1),1);
invalidDates = cell(1,size(raw,2));
anyInvalidDates = false(size(raw,1),1);
for col=[11,12]% Convert the contents of columns with dates to MATLAB datetimes using date format string.
    try
        dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[11,12]}, 'InputFormat', dateFormats{col==[11,12]}); %#ok<SAGROW>
    catch
        try
            % Handle dates surrounded by quotes
            dataArray{col} = cellfun(@(x) x(2:end-1), dataArray{col}, 'UniformOutput', false);
            dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[11,12]}, 'InputFormat', dateFormats{col==[11,12]}); %%#ok<SAGROW>
        catch
            dates{col} = repmat(datetime([NaN NaN NaN]), size(dataArray{col})); %#ok<SAGROW>
        end
    end
    
    dateFormatIndex = dateFormatIndex + 1;
    blankDates{col} = cellfun(@isempty, dataArray{col});
    anyBlankDates = blankDates{col} | anyBlankDates;
    invalidDates{col} = isnan(dates{col}.Hour) - blankDates{col};
    anyInvalidDates = invalidDates{col} | anyInvalidDates;
end
dates = dates(:,[11,12]);
blankDates = blankDates(:,[11,12]);
invalidDates = invalidDates(:,[11,12]);

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
maxdoas = table;
maxdoas.SpecNo = cell2mat(rawNumericColumns(:, 1));
maxdoas.Year = cell2mat(rawNumericColumns(:, 2));
maxdoas.Fractionalday = cell2mat(rawNumericColumns(:, 3));
maxdoas.Fractionaltime = cell2mat(rawNumericColumns(:, 4));
maxdoas.Scans = cell2mat(rawNumericColumns(:, 5));
maxdoas.Tint = cell2mat(rawNumericColumns(:, 6));
maxdoas.SZA = cell2mat(rawNumericColumns(:, 7));
maxdoas.SolarAzimuthAngle = cell2mat(rawNumericColumns(:, 8));
maxdoas.Elevviewingangle = cell2mat(rawNumericColumns(:, 9));
maxdoas.Azimviewingangle = cell2mat(rawNumericColumns(:, 10));
maxdoas.DateDDMMYYYY = dates{:, 1};
maxdoas.Timehhmmss = dates{:, 2};

maxdoas.Timehhmmss.Format='HH:mm:ss';

maxdoas.DateTime=maxdoas.DateDDMMYYYY+timeofday(maxdoas.Timehhmmss);
maxdoas.DateTime.Format='dd/MM/uuuu HH:mm:ss';

maxdoas.TotalExperimentTimesec = cell2mat(rawNumericColumns(:, 11));
maxdoas.NO2RMS = cell2mat(rawNumericColumns(:, 12));
maxdoas.NO2RefZm = cell2mat(rawNumericColumns(:, 13));
maxdoas.NO2processing_error = cell2mat(rawNumericColumns(:, 14));
maxdoas.NO2SlColo4 = cell2mat(rawNumericColumns(:, 15));
maxdoas.NO2SlErro4 = cell2mat(rawNumericColumns(:, 16));
maxdoas.NO2SlColRing = cell2mat(rawNumericColumns(:, 17));
maxdoas.NO2SlErrRing = cell2mat(rawNumericColumns(:, 18));
maxdoas.NO2SlColno2 = cell2mat(rawNumericColumns(:, 19));
maxdoas.NO2SlErrno2 = cell2mat(rawNumericColumns(:, 20));
maxdoas.NO2SlColo3 = cell2mat(rawNumericColumns(:, 21));
maxdoas.NO2SlErro3 = cell2mat(rawNumericColumns(:, 22));
maxdoas.NO2SlColoclo = cell2mat(rawNumericColumns(:, 23));
maxdoas.NO2SlErroclo = cell2mat(rawNumericColumns(:, 24));
maxdoas.NO2SlColbro = cell2mat(rawNumericColumns(:, 25));
maxdoas.NO2SlErrbro = cell2mat(rawNumericColumns(:, 26));
maxdoas.NO2SlColo3_243K = cell2mat(rawNumericColumns(:, 27));
maxdoas.NO2ShiftSpectrum = cell2mat(rawNumericColumns(:, 28));
maxdoas.NO2StretchSpectrum1 = cell2mat(rawNumericColumns(:, 29));
maxdoas.NO2StretchSpectrum2 = cell2mat(rawNumericColumns(:, 30));
maxdoas.BrORMS = cell2mat(rawNumericColumns(:, 31));
maxdoas.BrORefZm = cell2mat(rawNumericColumns(:, 32));
maxdoas.BrOprocessing_error = cell2mat(rawNumericColumns(:, 33));
maxdoas.BrOSlColo4 = cell2mat(rawNumericColumns(:, 34));
maxdoas.BrOSlErro4 = cell2mat(rawNumericColumns(:, 35));
maxdoas.BrOSlColRing = cell2mat(rawNumericColumns(:, 36));
maxdoas.BrOSlErrRing = cell2mat(rawNumericColumns(:, 37));
maxdoas.BrOSlColno2 = cell2mat(rawNumericColumns(:, 38));
maxdoas.BrOSlErrno2 = cell2mat(rawNumericColumns(:, 39));
maxdoas.BrOSlColo3 = cell2mat(rawNumericColumns(:, 40));
maxdoas.BrOSlErro3 = cell2mat(rawNumericColumns(:, 41));
maxdoas.BrOSlColoclo = cell2mat(rawNumericColumns(:, 42));
maxdoas.BrOSlErroclo = cell2mat(rawNumericColumns(:, 43));
maxdoas.BrOSlColbro = cell2mat(rawNumericColumns(:, 44));
maxdoas.BrOSlErrbro = cell2mat(rawNumericColumns(:, 45));
maxdoas.BrOSlColo3_243K = cell2mat(rawNumericColumns(:, 46));
maxdoas.BrOShiftSpectrum = cell2mat(rawNumericColumns(:, 47));
maxdoas.BrOStretchSpectrum1 = cell2mat(rawNumericColumns(:, 48));
maxdoas.BrOStretchSpectrum2 = cell2mat(rawNumericColumns(:, 49));
maxdoas.O4RMS = cell2mat(rawNumericColumns(:, 50));
maxdoas.O4RefZm = cell2mat(rawNumericColumns(:, 51));
maxdoas.O4processing_error = cell2mat(rawNumericColumns(:, 52));
maxdoas.O4SlColo4 = cell2mat(rawNumericColumns(:, 53));
maxdoas.O4SlErro4 = cell2mat(rawNumericColumns(:, 54));
maxdoas.O4SlColRing = cell2mat(rawNumericColumns(:, 55));
maxdoas.O4SlErrRing = cell2mat(rawNumericColumns(:, 56));
maxdoas.O4SlColno2 = cell2mat(rawNumericColumns(:, 57));
maxdoas.O4SlErrno2 = cell2mat(rawNumericColumns(:, 58));
maxdoas.O4SlColo3 = cell2mat(rawNumericColumns(:, 59));
maxdoas.O4SlErro3 = cell2mat(rawNumericColumns(:, 60));
maxdoas.O4SlColoclo = cell2mat(rawNumericColumns(:, 61));
maxdoas.O4SlErroclo = cell2mat(rawNumericColumns(:, 62));
maxdoas.O4SlColbro = cell2mat(rawNumericColumns(:, 63));
maxdoas.O4SlErrbro = cell2mat(rawNumericColumns(:, 64));
maxdoas.O4SlColo3_243K = cell2mat(rawNumericColumns(:, 65));
maxdoas.O4ShiftSpectrum = cell2mat(rawNumericColumns(:, 66));
maxdoas.O4StretchSpectrum1 = cell2mat(rawNumericColumns(:, 67));
maxdoas.O4StretchSpectrum2 = cell2mat(rawNumericColumns(:, 68));
maxdoas.Fluxes355 = cell2mat(rawNumericColumns(:, 69));
maxdoas.Fluxes360 = cell2mat(rawNumericColumns(:, 70));
maxdoas.Fluxes380 = cell2mat(rawNumericColumns(:, 71));
maxdoas.Fluxes385 = cell2mat(rawNumericColumns(:, 72));
maxdoas.Fluxes390 = cell2mat(rawNumericColumns(:, 73));
maxdoas.Fluxes405 = cell2mat(rawNumericColumns(:, 74));
maxdoas.Fluxes420 = cell2mat(rawNumericColumns(:, 75));
maxdoas.Fluxes425 = cell2mat(rawNumericColumns(:, 76));
maxdoas.Fluxes435 = cell2mat(rawNumericColumns(:, 77));
maxdoas.Fluxes440 = cell2mat(rawNumericColumns(:, 78));
maxdoas.Fluxes445 = cell2mat(rawNumericColumns(:, 79));
maxdoas.Fluxes450 = cell2mat(rawNumericColumns(:, 80));
maxdoas.Fluxes455 = cell2mat(rawNumericColumns(:, 81));
maxdoas.Fluxes460 = cell2mat(rawNumericColumns(:, 82));
maxdoas.Fluxes470 = cell2mat(rawNumericColumns(:, 83));
maxdoas.Fluxes490 = cell2mat(rawNumericColumns(:, 84));
maxdoas.Fluxes500 = cell2mat(rawNumericColumns(:, 85));
maxdoas.Fluxes532 = cell2mat(rawNumericColumns(:, 86));
maxdoas.Fluxes550 = cell2mat(rawNumericColumns(:, 87));

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% maxdoasd62151.DateDDMMYYYY=datenum(maxdoasd62151.DateDDMMYYYY);
% maxdoasd62151.Timehhmmss=datenum(maxdoasd62151.Timehhmmss);

% filter obviously bad data
maxdoas(maxdoas.BrOSlColbro>2e15,:)=[];

% save using same name
tmp=strsplit(filename,'.');
savename=[tmp{1} '.mat'];

save([path savename], 'maxdoas');


end
