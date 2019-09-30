function [  ] = replace_inst_label( )
%replace_inst_label: replace accidental PEARL-GBS header in UT-GBS files
%   
%   INPUT:  old csv files in ./old_csv/ directory
%   OUTPUT: files in directory defined by outfolder with az replaced


% Check if files are in the right place
if ~exist('old_csv', 'dir')
    error('Error: csv files must be placed in directory ./old_csv/');
end

list=dir('old_csv/');  % get info of files/folders in directory
isfile=~[list.isdir]; % determine index of files vs folders
filenames={list(isfile).name}; % create cell array of file names

% Set output directory (./ to create files in current folder)
outfolder='csv/';

% Create output folder if it doesn't exist, check if empty
if ~strcmp(outfolder,'./')
    if ~exist(outfolder, 'dir')
      mkdir(outfolder);
      fprintf('Created outout directory %s\n', outfolder);
      
    elseif size(dir(outfolder),1) > 2
        xx = input('Folder not empty, Press Enter to continue. ');
    end        
end


%%% Loop through files in old_csv/ and replace az values
for i=1:size(filenames,2)
    
    % Skip if not .csv file
    if isempty(strfind(filenames{i},'.csv'))
        fprintf('Skipping file %s\n', filenames{i});
        continue
    end
    
    % Set filenames
    oldpath = sprintf('old_csv/%s', filenames{i});
    newpath = strcat(outfolder, filenames{i});

    % Open files
    fid = fopen(oldpath,'r');
    newfile = fopen(newpath,'w');

    % Replace az value (assuming no other parameter is set to 180)
    old_content = fread(fid,'*char')';
    new_content = strrep(old_content, 'PEARL-GBS,', 'UT-GBS,');
    fprintf(newfile,'%s',new_content);

    fclose(fid);
    fclose(newfile);
    
    fprintf('Processed %s\n', filenames{i});

end

end