function data = read_prof_dscd(filename)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   read_prof_dscd = read_prof_dscd(FILENAME) Reads dSCD info for the given
%   day from HEIPro output
%
%   dSCDs are optical depth! Divide by cross-section at given wavelength
%
% Auto-generated by MATLAB on 2018/11/13 14:33:26
% modified by Kristof Bognar

%% Initialize variables.
delimiter = ' ';
startRow = 2;
endRow = inf;

%% Format string for each line of text:
%   column1: datetimes (%{dd/MM/yyyy}D)
%	column2: datetimes (%{HH:mm:ss}D)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%{dd/MM/yyyy}D%{HH:mm:ss}D%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% merge date and time

tmp=dataArray{1};
tmp2=dataArray{2};

tmp.Hour=hour(tmp2);
tmp.Minute=minute(tmp2);
tmp.Second=second(tmp2);

tmp.Format='dd/MM/yyyy HH:mm:ss';

tmp=table(tmp, 'VariableNames', {'date_time'});

%% Create output variable
data = table(dataArray{3:end-1}, 'VariableNames', {'SZA','elev','rel_azim','wl','meas','err_meas','retr'});

data = [tmp,data];


% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% meas20180324.date=datenum(meas20180324.date);
% meas20180324.time=datenum(meas20180324.time);

