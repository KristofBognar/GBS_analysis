function data = reformat_qdoas_RD(instr,batch)

year=datetime(now,'convertfrom','datenum').Year;
year=num2str(year);

batch_str=num2str(batch);

switch instr
    case 'u'
        filename=['/home/kristof/work/GBS/UT-GBS/' year '/QDOAS_output/UT-GBS_'...
                  year '_' batch_str '.ASC'];
        
        data = import_ut(filename);
        
        save(['/home/kristof/work/GBS/QDOAS_results/NDACC_RD_tables/UT-GBS_'...
              year '_' batch_str '.mat'],'data')
        
    case 'p'
        filename=['/home/kristof/work/GBS/PEARL-GBS/' year '/QDOAS_output/PEARL-GBS_'...
                  year '_' batch_str '.ASC'];
              
        data = import_p(filename);
              
        save(['/home/kristof/work/GBS/QDOAS_results/NDACC_RD_tables/PEARL-GBS_'...
              year '_' batch_str '.mat'],'data')
          
end
end


function data=import_ut(filename)
    %% Initialize variables.
    delimiter = '\t';
    startRow = 3;
    endRow = inf;

    %% Read columns of data as strings:
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

    %% Open the text file.
    fileID = fopen(filename,'r');

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

    for col=[1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246]
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
    rawNumericColumns = raw(:, [1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246]);

    %% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

    %% Create output variable
    data = table;
    data.SpecNo = cell2mat(rawNumericColumns(:, 1));
    data.Year = cell2mat(rawNumericColumns(:, 2));
    data.Fractionalday = cell2mat(rawNumericColumns(:, 3));
    data.Fractionaltime = cell2mat(rawNumericColumns(:, 4));
    data.Scans = cell2mat(rawNumericColumns(:, 5));
    data.Tint = cell2mat(rawNumericColumns(:, 6));
    data.SZA = cell2mat(rawNumericColumns(:, 7));
    data.SolarAzimuthAngle = cell2mat(rawNumericColumns(:, 8));
    data.Elevviewingangle = cell2mat(rawNumericColumns(:, 9));
    data.Azimviewingangle = cell2mat(rawNumericColumns(:, 10));
    data.DateDDMMYYYY = dates{:, 1};
    data.Timehhmmss = dates{:, 2};
    data.TotalExperimentTimesec = cell2mat(rawNumericColumns(:, 11));
    data.O4_VIS_293_a203RMS = cell2mat(rawNumericColumns(:, 12));
    data.O4_VIS_293_a203RefZm = cell2mat(rawNumericColumns(:, 13));
    data.O4_VIS_293_a203processing_error = cell2mat(rawNumericColumns(:, 14));
    data.O4_VIS_293_a203SlColno2 = cell2mat(rawNumericColumns(:, 15));
    data.O4_VIS_293_a203SlErrno2 = cell2mat(rawNumericColumns(:, 16));
    data.O4_VIS_293_a203SlColno2a = cell2mat(rawNumericColumns(:, 17));
    data.O4_VIS_293_a203SlErrno2a = cell2mat(rawNumericColumns(:, 18));
    data.O4_VIS_293_a203SlColo3 = cell2mat(rawNumericColumns(:, 19));
    data.O4_VIS_293_a203SlErro3 = cell2mat(rawNumericColumns(:, 20));
    data.O4_VIS_293_a203SlColo4 = cell2mat(rawNumericColumns(:, 21));
    data.O4_VIS_293_a203SlErro4 = cell2mat(rawNumericColumns(:, 22));
    data.O4_VIS_293_a203SlColo4a = cell2mat(rawNumericColumns(:, 23));
    data.O4_VIS_293_a203SlErro4a = cell2mat(rawNumericColumns(:, 24));
    data.O4_VIS_293_a203SlColh2o = cell2mat(rawNumericColumns(:, 25));
    data.O4_VIS_293_a203SlErrh2o = cell2mat(rawNumericColumns(:, 26));
    data.O4_VIS_293_a203SlColRing = cell2mat(rawNumericColumns(:, 27));
    data.O4_VIS_293_a203SlErrRing = cell2mat(rawNumericColumns(:, 28));
    data.O4_VIS_293_a203SlColx0 = cell2mat(rawNumericColumns(:, 29));
    data.O4_VIS_293_a203SlErrx0 = cell2mat(rawNumericColumns(:, 30));
    data.O4_VIS_293_a203SlColx1 = cell2mat(rawNumericColumns(:, 31));
    data.O4_VIS_293_a203SlErrx1 = cell2mat(rawNumericColumns(:, 32));
    data.O4_VIS_293_a203SlColx2 = cell2mat(rawNumericColumns(:, 33));
    data.O4_VIS_293_a203SlErrx2 = cell2mat(rawNumericColumns(:, 34));
    data.O4_VIS_293_a203SlColx3 = cell2mat(rawNumericColumns(:, 35));
    data.O4_VIS_293_a203SlErrx3 = cell2mat(rawNumericColumns(:, 36));
    data.O4_VIS_293_a203SlColx4 = cell2mat(rawNumericColumns(:, 37));
    data.O4_VIS_293_a203SlErrx4 = cell2mat(rawNumericColumns(:, 38));
    data.O4_VIS_293_a203SlColx5 = cell2mat(rawNumericColumns(:, 39));
    data.O4_VIS_293_a203SlErrx5 = cell2mat(rawNumericColumns(:, 40));
    data.O4_VIS_293_a203OffsetConstant = cell2mat(rawNumericColumns(:, 41));
    data.O4_VIS_293_a203ErrOffsetConstant = cell2mat(rawNumericColumns(:, 42));
    data.O4_VIS_293_a203ShiftSpectrum = cell2mat(rawNumericColumns(:, 43));
    data.O4_VIS_293_a203StretchSpectrum1 = cell2mat(rawNumericColumns(:, 44));
    data.O4_VIS_293_a203StretchSpectrum2 = cell2mat(rawNumericColumns(:, 45));
    data.O4_VIS_203_a293RMS = cell2mat(rawNumericColumns(:, 46));
    data.O4_VIS_203_a293RefZm = cell2mat(rawNumericColumns(:, 47));
    data.O4_VIS_203_a293processing_error = cell2mat(rawNumericColumns(:, 48));
    data.O4_VIS_203_a293SlColno2 = cell2mat(rawNumericColumns(:, 49));
    data.O4_VIS_203_a293SlErrno2 = cell2mat(rawNumericColumns(:, 50));
    data.O4_VIS_203_a293SlColno2a = cell2mat(rawNumericColumns(:, 51));
    data.O4_VIS_203_a293SlErrno2a = cell2mat(rawNumericColumns(:, 52));
    data.O4_VIS_203_a293SlColo3 = cell2mat(rawNumericColumns(:, 53));
    data.O4_VIS_203_a293SlErro3 = cell2mat(rawNumericColumns(:, 54));
    data.O4_VIS_203_a293SlColo4 = cell2mat(rawNumericColumns(:, 55));
    data.O4_VIS_203_a293SlErro4 = cell2mat(rawNumericColumns(:, 56));
    data.O4_VIS_203_a293SlColo4a = cell2mat(rawNumericColumns(:, 57));
    data.O4_VIS_203_a293SlErro4a = cell2mat(rawNumericColumns(:, 58));
    data.O4_VIS_203_a293SlColh2o = cell2mat(rawNumericColumns(:, 59));
    data.O4_VIS_203_a293SlErrh2o = cell2mat(rawNumericColumns(:, 60));
    data.O4_VIS_203_a293SlColRing = cell2mat(rawNumericColumns(:, 61));
    data.O4_VIS_203_a293SlErrRing = cell2mat(rawNumericColumns(:, 62));
    data.O4_VIS_203_a293SlColx0 = cell2mat(rawNumericColumns(:, 63));
    data.O4_VIS_203_a293SlErrx0 = cell2mat(rawNumericColumns(:, 64));
    data.O4_VIS_203_a293SlColx1 = cell2mat(rawNumericColumns(:, 65));
    data.O4_VIS_203_a293SlErrx1 = cell2mat(rawNumericColumns(:, 66));
    data.O4_VIS_203_a293SlColx2 = cell2mat(rawNumericColumns(:, 67));
    data.O4_VIS_203_a293SlErrx2 = cell2mat(rawNumericColumns(:, 68));
    data.O4_VIS_203_a293SlColx3 = cell2mat(rawNumericColumns(:, 69));
    data.O4_VIS_203_a293SlErrx3 = cell2mat(rawNumericColumns(:, 70));
    data.O4_VIS_203_a293SlColx4 = cell2mat(rawNumericColumns(:, 71));
    data.O4_VIS_203_a293SlErrx4 = cell2mat(rawNumericColumns(:, 72));
    data.O4_VIS_203_a293SlColx5 = cell2mat(rawNumericColumns(:, 73));
    data.O4_VIS_203_a293SlErrx5 = cell2mat(rawNumericColumns(:, 74));
    data.O4_VIS_203_a293OffsetConstant = cell2mat(rawNumericColumns(:, 75));
    data.O4_VIS_203_a293ErrOffsetConstant = cell2mat(rawNumericColumns(:, 76));
    data.O4_VIS_203_a293ShiftSpectrum = cell2mat(rawNumericColumns(:, 77));
    data.O4_VIS_203_a293StretchSpectrum1 = cell2mat(rawNumericColumns(:, 78));
    data.O4_VIS_203_a293StretchSpectrum2 = cell2mat(rawNumericColumns(:, 79));
    data.O4_VIS_203RMS = cell2mat(rawNumericColumns(:, 80));
    data.O4_VIS_203RefZm = cell2mat(rawNumericColumns(:, 81));
    data.O4_VIS_203processing_error = cell2mat(rawNumericColumns(:, 82));
    data.O4_VIS_203SlColno2 = cell2mat(rawNumericColumns(:, 83));
    data.O4_VIS_203SlErrno2 = cell2mat(rawNumericColumns(:, 84));
    data.O4_VIS_203SlColno2a = cell2mat(rawNumericColumns(:, 85));
    data.O4_VIS_203SlErrno2a = cell2mat(rawNumericColumns(:, 86));
    data.O4_VIS_203SlColo3 = cell2mat(rawNumericColumns(:, 87));
    data.O4_VIS_203SlErro3 = cell2mat(rawNumericColumns(:, 88));
    data.O4_VIS_203SlColo4 = cell2mat(rawNumericColumns(:, 89));
    data.O4_VIS_203SlErro4 = cell2mat(rawNumericColumns(:, 90));
    data.O4_VIS_203SlColh2o = cell2mat(rawNumericColumns(:, 91));
    data.O4_VIS_203SlErrh2o = cell2mat(rawNumericColumns(:, 92));
    data.O4_VIS_203SlColRing = cell2mat(rawNumericColumns(:, 93));
    data.O4_VIS_203SlErrRing = cell2mat(rawNumericColumns(:, 94));
    data.O4_VIS_203SlColx0 = cell2mat(rawNumericColumns(:, 95));
    data.O4_VIS_203SlErrx0 = cell2mat(rawNumericColumns(:, 96));
    data.O4_VIS_203SlColx1 = cell2mat(rawNumericColumns(:, 97));
    data.O4_VIS_203SlErrx1 = cell2mat(rawNumericColumns(:, 98));
    data.O4_VIS_203SlColx2 = cell2mat(rawNumericColumns(:, 99));
    data.O4_VIS_203SlErrx2 = cell2mat(rawNumericColumns(:, 100));
    data.O4_VIS_203SlColx3 = cell2mat(rawNumericColumns(:, 101));
    data.O4_VIS_203SlErrx3 = cell2mat(rawNumericColumns(:, 102));
    data.O4_VIS_203SlColx4 = cell2mat(rawNumericColumns(:, 103));
    data.O4_VIS_203SlErrx4 = cell2mat(rawNumericColumns(:, 104));
    data.O4_VIS_203SlColx5 = cell2mat(rawNumericColumns(:, 105));
    data.O4_VIS_203SlErrx5 = cell2mat(rawNumericColumns(:, 106));
    data.O4_VIS_203OffsetConstant = cell2mat(rawNumericColumns(:, 107));
    data.O4_VIS_203ErrOffsetConstant = cell2mat(rawNumericColumns(:, 108));
    data.O4_VIS_203ShiftSpectrum = cell2mat(rawNumericColumns(:, 109));
    data.O4_VIS_203StretchSpectrum1 = cell2mat(rawNumericColumns(:, 110));
    data.O4_VIS_203StretchSpectrum2 = cell2mat(rawNumericColumns(:, 111));
    data.O4_VIS_293RMS = cell2mat(rawNumericColumns(:, 112));
    data.O4_VIS_293RefZm = cell2mat(rawNumericColumns(:, 113));
    data.O4_VIS_293processing_error = cell2mat(rawNumericColumns(:, 114));
    data.O4_VIS_293SlColno2 = cell2mat(rawNumericColumns(:, 115));
    data.O4_VIS_293SlErrno2 = cell2mat(rawNumericColumns(:, 116));
    data.O4_VIS_293SlColno2a = cell2mat(rawNumericColumns(:, 117));
    data.O4_VIS_293SlErrno2a = cell2mat(rawNumericColumns(:, 118));
    data.O4_VIS_293SlColo3 = cell2mat(rawNumericColumns(:, 119));
    data.O4_VIS_293SlErro3 = cell2mat(rawNumericColumns(:, 120));
    data.O4_VIS_293SlColo4 = cell2mat(rawNumericColumns(:, 121));
    data.O4_VIS_293SlErro4 = cell2mat(rawNumericColumns(:, 122));
    data.O4_VIS_293SlColh2o = cell2mat(rawNumericColumns(:, 123));
    data.O4_VIS_293SlErrh2o = cell2mat(rawNumericColumns(:, 124));
    data.O4_VIS_293SlColRing = cell2mat(rawNumericColumns(:, 125));
    data.O4_VIS_293SlErrRing = cell2mat(rawNumericColumns(:, 126));
    data.O4_VIS_293SlColx0 = cell2mat(rawNumericColumns(:, 127));
    data.O4_VIS_293SlErrx0 = cell2mat(rawNumericColumns(:, 128));
    data.O4_VIS_293SlColx1 = cell2mat(rawNumericColumns(:, 129));
    data.O4_VIS_293SlErrx1 = cell2mat(rawNumericColumns(:, 130));
    data.O4_VIS_293SlColx2 = cell2mat(rawNumericColumns(:, 131));
    data.O4_VIS_293SlErrx2 = cell2mat(rawNumericColumns(:, 132));
    data.O4_VIS_293SlColx3 = cell2mat(rawNumericColumns(:, 133));
    data.O4_VIS_293SlErrx3 = cell2mat(rawNumericColumns(:, 134));
    data.O4_VIS_293SlColx4 = cell2mat(rawNumericColumns(:, 135));
    data.O4_VIS_293SlErrx4 = cell2mat(rawNumericColumns(:, 136));
    data.O4_VIS_293SlColx5 = cell2mat(rawNumericColumns(:, 137));
    data.O4_VIS_293SlErrx5 = cell2mat(rawNumericColumns(:, 138));
    data.O4_VIS_293OffsetConstant = cell2mat(rawNumericColumns(:, 139));
    data.O4_VIS_293ErrOffsetConstant = cell2mat(rawNumericColumns(:, 140));
    data.O4_VIS_293ShiftSpectrum = cell2mat(rawNumericColumns(:, 141));
    data.O4_VIS_293StretchSpectrum1 = cell2mat(rawNumericColumns(:, 142));
    data.O4_VIS_293StretchSpectrum2 = cell2mat(rawNumericColumns(:, 143));
    data.O3_noXRMS = cell2mat(rawNumericColumns(:, 144));
    data.O3_noXRefZm = cell2mat(rawNumericColumns(:, 145));
    data.O3_noXprocessing_error = cell2mat(rawNumericColumns(:, 146));
    data.O3_noXSlColh2o = cell2mat(rawNumericColumns(:, 147));
    data.O3_noXSlErrh2o = cell2mat(rawNumericColumns(:, 148));
    data.O3_noXSlColo4 = cell2mat(rawNumericColumns(:, 149));
    data.O3_noXSlErro4 = cell2mat(rawNumericColumns(:, 150));
    data.O3_noXSlColRing = cell2mat(rawNumericColumns(:, 151));
    data.O3_noXSlErrRing = cell2mat(rawNumericColumns(:, 152));
    data.O3_noXSlColno2 = cell2mat(rawNumericColumns(:, 153));
    data.O3_noXSlErrno2 = cell2mat(rawNumericColumns(:, 154));
    data.O3_noXSlColo3 = cell2mat(rawNumericColumns(:, 155));
    data.O3_noXSlErro3 = cell2mat(rawNumericColumns(:, 156));
    data.O3_noXSlColx0 = cell2mat(rawNumericColumns(:, 157));
    data.O3_noXSlErrx0 = cell2mat(rawNumericColumns(:, 158));
    data.O3_noXSlColx1 = cell2mat(rawNumericColumns(:, 159));
    data.O3_noXSlErrx1 = cell2mat(rawNumericColumns(:, 160));
    data.O3_noXSlColx2 = cell2mat(rawNumericColumns(:, 161));
    data.O3_noXSlErrx2 = cell2mat(rawNumericColumns(:, 162));
    data.O3_noXSlColx3 = cell2mat(rawNumericColumns(:, 163));
    data.O3_noXSlErrx3 = cell2mat(rawNumericColumns(:, 164));
    data.O3_noXShiftSpectrum = cell2mat(rawNumericColumns(:, 165));
    data.O3_noXStretchSpectrum1 = cell2mat(rawNumericColumns(:, 166));
    data.O3_noXStretchSpectrum2 = cell2mat(rawNumericColumns(:, 167));
    data.O3RMS = cell2mat(rawNumericColumns(:, 168));
    data.O3RefZm = cell2mat(rawNumericColumns(:, 169));
    data.O3processing_error = cell2mat(rawNumericColumns(:, 170));
    data.O3SlColh2o = cell2mat(rawNumericColumns(:, 171));
    data.O3SlErrh2o = cell2mat(rawNumericColumns(:, 172));
    data.O3SlColo4 = cell2mat(rawNumericColumns(:, 173));
    data.O3SlErro4 = cell2mat(rawNumericColumns(:, 174));
    data.O3SlColRing = cell2mat(rawNumericColumns(:, 175));
    data.O3SlErrRing = cell2mat(rawNumericColumns(:, 176));
    data.O3SlColno2 = cell2mat(rawNumericColumns(:, 177));
    data.O3SlErrno2 = cell2mat(rawNumericColumns(:, 178));
    data.O3SlColo3 = cell2mat(rawNumericColumns(:, 179));
    data.O3SlErro3 = cell2mat(rawNumericColumns(:, 180));
    data.O3SlColX = cell2mat(rawNumericColumns(:, 181));
    data.O3SlErrX = cell2mat(rawNumericColumns(:, 182));
    data.O3SlColx0 = cell2mat(rawNumericColumns(:, 183));
    data.O3SlErrx0 = cell2mat(rawNumericColumns(:, 184));
    data.O3SlColx1 = cell2mat(rawNumericColumns(:, 185));
    data.O3SlErrx1 = cell2mat(rawNumericColumns(:, 186));
    data.O3SlColx2 = cell2mat(rawNumericColumns(:, 187));
    data.O3SlErrx2 = cell2mat(rawNumericColumns(:, 188));
    data.O3SlColx3 = cell2mat(rawNumericColumns(:, 189));
    data.O3SlErrx3 = cell2mat(rawNumericColumns(:, 190));
    data.O3ShiftSpectrum = cell2mat(rawNumericColumns(:, 191));
    data.O3StretchSpectrum1 = cell2mat(rawNumericColumns(:, 192));
    data.O3StretchSpectrum2 = cell2mat(rawNumericColumns(:, 193));
    data.NO2RMS = cell2mat(rawNumericColumns(:, 194));
    data.NO2RefZm = cell2mat(rawNumericColumns(:, 195));
    data.NO2processing_error = cell2mat(rawNumericColumns(:, 196));
    data.NO2SlColh2o = cell2mat(rawNumericColumns(:, 197));
    data.NO2SlErrh2o = cell2mat(rawNumericColumns(:, 198));
    data.NO2SlColo4 = cell2mat(rawNumericColumns(:, 199));
    data.NO2SlErro4 = cell2mat(rawNumericColumns(:, 200));
    data.NO2SlColRing = cell2mat(rawNumericColumns(:, 201));
    data.NO2SlErrRing = cell2mat(rawNumericColumns(:, 202));
    data.NO2SlColno2 = cell2mat(rawNumericColumns(:, 203));
    data.NO2SlErrno2 = cell2mat(rawNumericColumns(:, 204));
    data.NO2SlColo3 = cell2mat(rawNumericColumns(:, 205));
    data.NO2SlErro3 = cell2mat(rawNumericColumns(:, 206));
    data.NO2ShiftSpectrum = cell2mat(rawNumericColumns(:, 207));
    data.NO2StretchSpectrum1 = cell2mat(rawNumericColumns(:, 208));
    data.NO2StretchSpectrum2 = cell2mat(rawNumericColumns(:, 209));
    data.NO2_425450RMS = cell2mat(rawNumericColumns(:, 210));
    data.NO2_425450RefZm = cell2mat(rawNumericColumns(:, 211));
    data.NO2_425450processing_error = cell2mat(rawNumericColumns(:, 212));
    data.NO2_425450SlColh2o = cell2mat(rawNumericColumns(:, 213));
    data.NO2_425450SlErrh2o = cell2mat(rawNumericColumns(:, 214));
    data.NO2_425450SlColo4 = cell2mat(rawNumericColumns(:, 215));
    data.NO2_425450SlErro4 = cell2mat(rawNumericColumns(:, 216));
    data.NO2_425450SlColRing = cell2mat(rawNumericColumns(:, 217));
    data.NO2_425450SlErrRing = cell2mat(rawNumericColumns(:, 218));
    data.NO2_425450SlColno2 = cell2mat(rawNumericColumns(:, 219));
    data.NO2_425450SlErrno2 = cell2mat(rawNumericColumns(:, 220));
    data.NO2_425450SlColo3 = cell2mat(rawNumericColumns(:, 221));
    data.NO2_425450SlErro3 = cell2mat(rawNumericColumns(:, 222));
    data.NO2_425450ShiftSpectrum = cell2mat(rawNumericColumns(:, 223));
    data.NO2_425450StretchSpectrum1 = cell2mat(rawNumericColumns(:, 224));
    data.NO2_425450StretchSpectrum2 = cell2mat(rawNumericColumns(:, 225));
    data.Fluxes355 = cell2mat(rawNumericColumns(:, 226));
    data.Fluxes360 = cell2mat(rawNumericColumns(:, 227));
    data.Fluxes380 = cell2mat(rawNumericColumns(:, 228));
    data.Fluxes385 = cell2mat(rawNumericColumns(:, 229));
    data.Fluxes390 = cell2mat(rawNumericColumns(:, 230));
    data.Fluxes405 = cell2mat(rawNumericColumns(:, 231));
    data.Fluxes420 = cell2mat(rawNumericColumns(:, 232));
    data.Fluxes425 = cell2mat(rawNumericColumns(:, 233));
    data.Fluxes435 = cell2mat(rawNumericColumns(:, 234));
    data.Fluxes440 = cell2mat(rawNumericColumns(:, 235));
    data.Fluxes445 = cell2mat(rawNumericColumns(:, 236));
    data.Fluxes450 = cell2mat(rawNumericColumns(:, 237));
    data.Fluxes455 = cell2mat(rawNumericColumns(:, 238));
    data.Fluxes460 = cell2mat(rawNumericColumns(:, 239));
    data.Fluxes470 = cell2mat(rawNumericColumns(:, 240));
    data.Fluxes490 = cell2mat(rawNumericColumns(:, 241));
    data.Fluxes500 = cell2mat(rawNumericColumns(:, 242));
    data.Fluxes532 = cell2mat(rawNumericColumns(:, 243));
    data.Fluxes550 = cell2mat(rawNumericColumns(:, 244));

    % For code requiring serial dates (datenum) instead of datetime, uncomment
    % the following line(s) below to return the imported dates as datenum(s).

    % UTGBS20193.DateDDMMYYYY=datenum(UTGBS20193.DateDDMMYYYY);
    % UTGBS20193.Timehhmmss=datenum(UTGBS20193.Timehhmmss);
    
end

function import_p(filename)
    %% Initialize variables.
    delimiter = '\t';
    startRow = 3;
    endRow = inf;

    %% Read columns of data as strings:
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

    %% Open the text file.
    fileID = fopen(filename,'r');

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

    for col=[1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78]
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
    rawNumericColumns = raw(:, [1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78]);

    %% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

    %% Create output variable
    PEARLGBS20192 = table;
    PEARLGBS20192.SpecNo = cell2mat(rawNumericColumns(:, 1));
    PEARLGBS20192.Year = cell2mat(rawNumericColumns(:, 2));
    PEARLGBS20192.Fractionalday = cell2mat(rawNumericColumns(:, 3));
    PEARLGBS20192.Fractionaltime = cell2mat(rawNumericColumns(:, 4));
    PEARLGBS20192.Scans = cell2mat(rawNumericColumns(:, 5));
    PEARLGBS20192.Tint = cell2mat(rawNumericColumns(:, 6));
    PEARLGBS20192.SZA = cell2mat(rawNumericColumns(:, 7));
    PEARLGBS20192.SolarAzimuthAngle = cell2mat(rawNumericColumns(:, 8));
    PEARLGBS20192.Elevviewingangle = cell2mat(rawNumericColumns(:, 9));
    PEARLGBS20192.Azimviewingangle = cell2mat(rawNumericColumns(:, 10));
    PEARLGBS20192.DateDDMMYYYY = dates{:, 1};
    PEARLGBS20192.Timehhmmss = dates{:, 2};
    PEARLGBS20192.TotalExperimentTimesec = cell2mat(rawNumericColumns(:, 11));
    PEARLGBS20192.O3RMS = cell2mat(rawNumericColumns(:, 12));
    PEARLGBS20192.O3RefZm = cell2mat(rawNumericColumns(:, 13));
    PEARLGBS20192.O3processing_error = cell2mat(rawNumericColumns(:, 14));
    PEARLGBS20192.O3SlColRing = cell2mat(rawNumericColumns(:, 15));
    PEARLGBS20192.O3SlErrRing = cell2mat(rawNumericColumns(:, 16));
    PEARLGBS20192.O3SlColno2 = cell2mat(rawNumericColumns(:, 17));
    PEARLGBS20192.O3SlErrno2 = cell2mat(rawNumericColumns(:, 18));
    PEARLGBS20192.O3SlColo3 = cell2mat(rawNumericColumns(:, 19));
    PEARLGBS20192.O3SlErro3 = cell2mat(rawNumericColumns(:, 20));
    PEARLGBS20192.O3SlColo3a = cell2mat(rawNumericColumns(:, 21));
    PEARLGBS20192.O3SlErro3a = cell2mat(rawNumericColumns(:, 22));
    PEARLGBS20192.O3SlColo3p1 = cell2mat(rawNumericColumns(:, 23));
    PEARLGBS20192.O3SlErro3p1 = cell2mat(rawNumericColumns(:, 24));
    PEARLGBS20192.O3SlColo3p2 = cell2mat(rawNumericColumns(:, 25));
    PEARLGBS20192.O3SlErro3p2 = cell2mat(rawNumericColumns(:, 26));
    PEARLGBS20192.O3SlColhcho = cell2mat(rawNumericColumns(:, 27));
    PEARLGBS20192.O3SlErrhcho = cell2mat(rawNumericColumns(:, 28));
    PEARLGBS20192.O3SlColx0 = cell2mat(rawNumericColumns(:, 29));
    PEARLGBS20192.O3SlErrx0 = cell2mat(rawNumericColumns(:, 30));
    PEARLGBS20192.O3SlColx1 = cell2mat(rawNumericColumns(:, 31));
    PEARLGBS20192.O3SlErrx1 = cell2mat(rawNumericColumns(:, 32));
    PEARLGBS20192.O3SlColx2 = cell2mat(rawNumericColumns(:, 33));
    PEARLGBS20192.O3SlErrx2 = cell2mat(rawNumericColumns(:, 34));
    PEARLGBS20192.O3SlColx3 = cell2mat(rawNumericColumns(:, 35));
    PEARLGBS20192.O3SlErrx3 = cell2mat(rawNumericColumns(:, 36));
    PEARLGBS20192.O3ShiftSpectrum = cell2mat(rawNumericColumns(:, 37));
    PEARLGBS20192.O3StretchSpectrum1 = cell2mat(rawNumericColumns(:, 38));
    PEARLGBS20192.O3StretchSpectrum2 = cell2mat(rawNumericColumns(:, 39));
    PEARLGBS20192.NO2RMS = cell2mat(rawNumericColumns(:, 40));
    PEARLGBS20192.NO2RefZm = cell2mat(rawNumericColumns(:, 41));
    PEARLGBS20192.NO2processing_error = cell2mat(rawNumericColumns(:, 42));
    PEARLGBS20192.NO2SlColo4 = cell2mat(rawNumericColumns(:, 43));
    PEARLGBS20192.NO2SlErro4 = cell2mat(rawNumericColumns(:, 44));
    PEARLGBS20192.NO2SlColRing = cell2mat(rawNumericColumns(:, 45));
    PEARLGBS20192.NO2SlErrRing = cell2mat(rawNumericColumns(:, 46));
    PEARLGBS20192.NO2SlColno2 = cell2mat(rawNumericColumns(:, 47));
    PEARLGBS20192.NO2SlErrno2 = cell2mat(rawNumericColumns(:, 48));
    PEARLGBS20192.NO2SlColo3 = cell2mat(rawNumericColumns(:, 49));
    PEARLGBS20192.NO2SlErro3 = cell2mat(rawNumericColumns(:, 50));
    PEARLGBS20192.NO2SlColoclo = cell2mat(rawNumericColumns(:, 51));
    PEARLGBS20192.NO2SlErroclo = cell2mat(rawNumericColumns(:, 52));
    PEARLGBS20192.NO2SlColbro = cell2mat(rawNumericColumns(:, 53));
    PEARLGBS20192.NO2SlErrbro = cell2mat(rawNumericColumns(:, 54));
    PEARLGBS20192.NO2ShiftSpectrum = cell2mat(rawNumericColumns(:, 55));
    PEARLGBS20192.NO2StretchSpectrum1 = cell2mat(rawNumericColumns(:, 56));
    PEARLGBS20192.NO2StretchSpectrum2 = cell2mat(rawNumericColumns(:, 57));
    PEARLGBS20192.Fluxes355 = cell2mat(rawNumericColumns(:, 58));
    PEARLGBS20192.Fluxes360 = cell2mat(rawNumericColumns(:, 59));
    PEARLGBS20192.Fluxes380 = cell2mat(rawNumericColumns(:, 60));
    PEARLGBS20192.Fluxes385 = cell2mat(rawNumericColumns(:, 61));
    PEARLGBS20192.Fluxes390 = cell2mat(rawNumericColumns(:, 62));
    PEARLGBS20192.Fluxes405 = cell2mat(rawNumericColumns(:, 63));
    PEARLGBS20192.Fluxes420 = cell2mat(rawNumericColumns(:, 64));
    PEARLGBS20192.Fluxes425 = cell2mat(rawNumericColumns(:, 65));
    PEARLGBS20192.Fluxes435 = cell2mat(rawNumericColumns(:, 66));
    PEARLGBS20192.Fluxes440 = cell2mat(rawNumericColumns(:, 67));
    PEARLGBS20192.Fluxes445 = cell2mat(rawNumericColumns(:, 68));
    PEARLGBS20192.Fluxes450 = cell2mat(rawNumericColumns(:, 69));
    PEARLGBS20192.Fluxes455 = cell2mat(rawNumericColumns(:, 70));
    PEARLGBS20192.Fluxes460 = cell2mat(rawNumericColumns(:, 71));
    PEARLGBS20192.Fluxes470 = cell2mat(rawNumericColumns(:, 72));
    PEARLGBS20192.Fluxes490 = cell2mat(rawNumericColumns(:, 73));
    PEARLGBS20192.Fluxes500 = cell2mat(rawNumericColumns(:, 74));
    PEARLGBS20192.Fluxes532 = cell2mat(rawNumericColumns(:, 75));
    PEARLGBS20192.Fluxes550 = cell2mat(rawNumericColumns(:, 76));

    % For code requiring serial dates (datenum) instead of datetime, uncomment
    % the following line(s) below to return the imported dates as datenum(s).

    % PEARLGBS20192.DateDDMMYYYY=datenum(PEARLGBS20192.DateDDMMYYYY);
    % PEARLGBS20192.Timehhmmss=datenum(PEARLGBS20192.Timehhmmss);

end
