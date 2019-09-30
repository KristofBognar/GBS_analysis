function read_qdoas_table(instr,readO3X,startyear,endyear,UV,merge_ut)
%read_qdoas_table(instr,readO3X,startyear,endyear,merge_ut)
% To read QDOAS output in table format and save yearly files
% File names and directories are hardcoded in the function
% File format (number and names of columns) is hardcoded!
%
% INPUT:    instr: 1 for UT-GBS, 2 for PEARL-GBS
%           readO3X: true or false, to read UT-GBS O3 data with X.xs
%           startyear,endyear: start and end years of iteration
%           UV (optional): boolean, reads UV data from both intruments if true
%           merge_ut (optional): merges UT-GBS datasets using rows common
%               to both tables -- O3X data must be read before using this option
% OUTPUT: saved yearly files 
%
% Created by Kristof Bognar, October 2017


% % % debug
% % % start and end of loop
% % startyear=1999;
% % endyear=2016;
% % 
% % % select UT-GBS (1) or PEARL-GBS (2)
% % instr=1;
% % 
% % % read files with O3 + X.cs fit for UT-GBS?
% % % if false, regular output (no X.xs) is read with all the windows (4*O4, O3, NO2)
% % readO3X=true;

%% input parameters

if nargin<6, merge_ut=false; end
if nargin<5, UV=false; end

O3X_str='';
UV_str='';
if instr==1
    if readO3X
        ftype=2;
        O3X_str='_O3X';
    else
        ftype=1;
    end
    
    if UV, ftype=4; O3X_str=''; UV_str='_UV'; end

elseif instr==2 
    ftype=3;
    if readO3X, warning('Only UT-GBS (instr 1) has separate O3X files'); end
    if UV, ftype=5; UV_str='_UV'; end
end



qdoas_dir='/home/kristof/work/GBS/QDOAS_results/';

% list of file names
if instr==1
    f_list={...
    'UT-GBS_1999_reanalysis_VIS_300',...
    'UT-GBS_2000_reanalysis_VIS',...
    'UT-GBS_2003_reanalysis_VIS_600+300',...
    'UT-GBS_2004_reanalysis_VIS_300',...
    'UT-GBS_2005_reanalysis_VIS',...
    'UT-GBS_2006_reanalysis_VIS',...
    'UT-GBS_2007_reanalysis_VIS',...
    'UT-GBS_2008_reanalysis_VIS',...
    'UT-GBS_2009_reanalysis_VIS',...
    'UT-GBS_2010_reanalysis_VIS',...
    'UT-GBS_2011_reanalysis_VIS',...
    'UT-GBS_2012_reanalysis_VIS',...
    'UT-GBS_2013_reanalysis_VIS',...
    'UT-GBS_2014_reanalysis_VIS',...
    'UT-GBS_2015_reanalysis_VIS',...
	'UT-GBS_2016_reanalysis_VIS',...
    'UT-GBS_2017_reanalysis_VIS'};

    if UV
        f_list={...
        'UT-GBS_2008_reanalysis_UV',...
        'UT-GBS_2009_reanalysis_UV',...
        'UT-GBS_2010_reanalysis_UV'};
    end
        
else
    f_list={...
    'PEARL-GBS_2006_reanalysis_VIS',...
    'PEARL-GBS_2007_reanalysis_VIS',...
    'PEARL-GBS_2008_reanalysis_VIS',...
    'PEARL-GBS_2009_reanalysis_VIS',...
    'PEARL-GBS_2010_reanalysis_VIS'};

    if UV
        f_list={...
        'PEARL-GBS_2007_reanalysis_UV',...
        'PEARL-GBS_2008_reanalysis_UV',...
        'PEARL-GBS_2009_reanalysis_UV',...
        'PEARL-GBS_2010_reanalysis_UV',...
        'PEARL-GBS_2011_reanalysis_UV',...
        'PEARL-GBS_2012_reanalysis_UV',...
        'PEARL-GBS_2013_reanalysis_UV',...
        'PEARL-GBS_2014_reanalysis_UV',...
        'PEARL-GBS_2015_reanalysis_UV',...
        'PEARL-GBS_2016_reanalysis_UV',...
        'PEARL-GBS_2017_reanalysis_UV'};
    end
        
end

% list of years 
years=NaN(length(f_list),1);
for i=1:length(f_list)
    tmp=strsplit(f_list{i},'_');
    years(i)=str2double(tmp{2});
end
% get instrument name
instr_str=tmp{1};
    
%% loop over desired years
n=0;
for i=max(min(years),startyear):min(max(years),endyear)
    
    % find file index
    ind=find(years==i);
    if isempty(ind), continue, end
    
    % Data format for UT-GBS changes starting from 2017
    if i>=2017 && instr==1 && ~UV
        ftype=3;
        O3X_str='';
    end

    % display progress info
    disp_str=['Reading ' num2str(i) ' ' instr_str O3X_str UV_str ' data'];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str);    
    
    % create full filename
    fname=[qdoas_dir f_list{ind} O3X_str '.ASC'];
    
    % read file
    data = read_qdoas_table__read_file(fname,ftype);

    % save yearly file
    savename=[qdoas_dir 'yearly_tables/' instr_str '_' num2str(i) O3X_str UV_str '.mat'];
    save(savename,'data');

    % merge UT-GBS datasets if option is selected
    if instr==1 && ~readO3X && merge_ut
        
        % save table
        o3=data;
        
        try % loading O3X file
            load([qdoas_dir 'yearly_tables/' instr_str '_' num2str(i) '_O3X.mat']);
            o3x=data;
        catch
            warning('No O3X data saved for this year, skipping merging');
            continue
        end
        
        if abs(length(o3.Year)-length(o3x.Year))>10
            size_diff=length(o3.Year)-length(o3x.Year);
            warning(['O3 - O3X size difference is ' num2str(size_diff) ' rows']);
        end
        
            % find common lines (should be near exact match)
        [~,io3,io3x]=intersect(o3.Fractionalday,o3x.Fractionalday);

        % clip tables to size
        o3=o3(io3,:);
        o3x=o3x(io3x,:);
        
        % merge tables
        data=[o3(:,1:185) o3x(:,14:39) o3(:,186:end)];
        
        % save merged table
        savename=[qdoas_dir 'yearly_tables/' instr_str '_' num2str(i) '_merged.mat'];
        save(savename,'data');
        
        
    end
    
end

fprintf('\n');
fprintf('Done\n');

