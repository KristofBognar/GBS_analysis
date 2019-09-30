% read specified column(s) from QDOAS reanalysis files and save them
% in map object with the year as the key


%% input parameters

% select UT-GBS (1) or PEARL-GBS (2)
instr=1;

%read visible or UV data?
read_UV=false;

% read files with O3 + X.cs fit for UT-GBS?
% if false, regular output (no X.xs) is read with all the windows (4*O4, O3, NO2)
readO3X=false;

if ~read_UV

    % create map object to hold yearly data
    if instr==1
        if readO3X
            utgbs_vis_o3x=containers.Map;
        else
            utgbs_vis=containers.Map;
        end
    elseif instr==2
        pgbs_vis=containers.Map;
    end
        
    % list of file names
    if instr==1 % UT-GBS
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
        'UT-GBS_2016_reanalysis_VIS'};

    elseif instr==2 % PEARL-GBS
        f_list={...
        'PEARL-GBS_2006_reanalysis_VIS',...
        'PEARL-GBS_2007_reanalysis_VIS',...
        'PEARL-GBS_2008_reanalysis_VIS',...
        'PEARL-GBS_2009_reanalysis_VIS',...
        'PEARL-GBS_2010_reanalysis_VIS'};
    end

    % get year from file names
    for i=1:size(f_list,2)
        tmp=strsplit(char(f_list{i}),'_');
        year{i}=tmp{2};
    end
    

    %% loop over the years
    tmp=1;
    for i=1:length(f_list)

        % create full filename
        if readO3X
            fname=[f_list{i} '_O3X.ASC'];
        else
            fname=[f_list{i} '.ASC'];
        end

        % progress info
        disp(fname);

        % read file header
        fid=fopen(fname);
        fgets(fid);
        line = fgets(fid);
        fclose(fid);

        header=strsplit(line,'\t');
        
        % make sure the header is the same as the previous file
        if i>1, tmp=length(header); end
        if tmp~=length(header) && i>1
            error('Error reading file: file format might have changed');
        end

        % read in full data file
        data=read_qdoas_map__read_file(fname,header);
        
        % assign data to map
        if instr==1
            if readO3X
                utgbs_vis_o3x(year{i})=data;
            else
                utgbs_vis(year{i})=data;
            end
        elseif instr==2
            pgbs_vis(year{i})=data;
        end
        

    end

    if instr==1
        if readO3X
            save map_UT-GBS_1999-2016_reanalysis_VIS_O3X.mat utgbs_vis_o3x
        else
            save map_UT-GBS_1999-2016_reanalysis_VIS.mat utgbs_vis
        end
    elseif instr==2
        save map_PEARL-GBS_2006-2010_reanalysis_VIS.mat pgbs_vis
    end
    
else
    
    % create map object to hold yearly data
    if instr==1
        utgbs_uv=containers.Map;
    elseif instr==2
        pgbs_uv=containers.Map;
    end
    
    % list of file names
    if instr==1
        f_list={...
                'UT-GBS_2008_reanalysis_UV.ASC'...
                'UT-GBS_2009_reanalysis_UV.ASC'...
                'UT-GBS_2010_reanalysis_UV.ASC'};
    elseif instr==2
        f_list={...
                'PEARL-GBS_2007_reanalysis_UV.ASC'...
                'PEARL-GBS_2008_reanalysis_UV.ASC'...
                'PEARL-GBS_2009_reanalysis_UV.ASC'...
                'PEARL-GBS_2010_reanalysis_UV.ASC'...
                'PEARL-GBS_2011_reanalysis_UV.ASC'...
                'PEARL-GBS_2012_reanalysis_UV.ASC'...
                'PEARL-GBS_2013_reanalysis_UV.ASC'...
                'PEARL-GBS_2014_reanalysis_UV.ASC'...
                'PEARL-GBS_2015_reanalysis_UV.ASC'...
                'PEARL-GBS_2016_reanalysis_UV.ASC'};
    end

    % get year from file names
    for i=1:size(f_list,2)
        tmp=strsplit(char(f_list{i}),'_');
        year{i}=tmp{2};
    end
    
    
    tmp=1;
    for i=1:length(f_list)

        fname=f_list{i};
        % progress info
        disp(fname);

        % read file header
        fid=fopen(fname);
        fgets(fid);
        line = fgets(fid);
        fclose(fid);

        if i>1, tmp=length(uv_header); end
        uv_header=strsplit(line,'\t');

        % make sure the header is the same as the previous file
        if tmp~=length(uv_header) && i>1
            error('Error reading file: file format might have changed');
        end

        % read in full data file
        data=read_qdoas_map__read_file(fname,uv_header);

        % assign data to map
        if instr==1
            utgbs_uv(year{i})=data;
        elseif instr==2
            pgbs_uv(year{i})=data;
        end


    end

    if instr==1
        save map_UT-GBS_2008-2010_reanalysis_UV.mat utgbs_uv
    elseif instr==2
        save map_PEARL-GBS_2007-2016_reanalysis_UV.mat pgbs_uv
    end

end


