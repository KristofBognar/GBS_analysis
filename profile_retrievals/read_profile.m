function profile = read_profile(filename,option)
%IMPORTFILE Import numeric data from a text file as a matrix.
% imports profile files from either aerosol or tracegas retrieval results
%
% INPUT:
%   filename: name of file to read
%   option: 'a' for aerosol files, 'tg' for tracegas files
%
% Auto-generated by MATLAB on 2017/02/01 14:19:28

%% Initialize variables.
startRow = 2;
endRow = inf;

%% Format string for each line of text:
% For more information, see the TEXTSCAN documentation.
if option=='a'
    formatSpec = '%14f%14f%14f%14f%14f%14f%f%[^\n\r]';
elseif option=='tg'
    formatSpec = '%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%f%[^\n\r]';
end
%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
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

%% Create output variable
profile = [dataArray{1:end-1}];