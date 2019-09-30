function [time_info, spec_info, sp] = CCD_tests(time_info, spec_info, sp)

    if nargin==0
    
        filename='/home/kristof/work/GBS/PEARL-GBS/2019/testing_Eureka/csv/01031919.csv';

        version=2;
        sat_lev = 60000;
        nbr_pix = 2048;
        date_format = 'dd/mm/yyyy';


        [time_info, spec_info, sp] = ...
        read_csv(filename, version, sat_lev, nbr_pix, date_format);   

        return
        
    end
    
    %% spectral noise test

    
    
    rms=[];
    specs=sp(8:62,:); % 100 ms, 1 spectrum each
    
    
    for i=1:size(specs,1)
        tmp=mean(specs(1:i,:),1);
        rms(i)=sqrt(mean(tmp.^2));
    end
    
    plot(rms), hold on
%     loglog([1:length(rms)].^(-0.5),'r--')
    
    
%     rms=[];
%     specs=sp(8:62,:); % 100 ms, 1 spectrum each
%     
%     for i=1:size(specs,1)
%         specs(i,:)=specs(i,:)/mean(specs(i,:));
%     end
%     
%     
%     tmp_mean=mean(specs,1);
%     
%     for i=1:size(specs,1)
%         tmp=mean(specs(1:i,:),1);
%         rms(i)=sqrt(mean((tmp./tmp_mean).^2));
%     end
%     
%     plot(rms), hold on
%     % plot([1:length(rms)].^(-0.5),'r--')
    
    
    
    


end

function [time_info, spec_info, sp, is_sat, version, nbr_cols] = ....
    read_csv(filename, version, sat_lev, nbr_pix, date_format)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Taken from Preprocs code, updated to return relevant info only, in table
% format
%
% This function reads all the spectra in a given csv file into lines that
% will correspond with the lines to be printed in the L1 and L2 files.
% This calculates the elevation angles for each given spectrum.
%
% Input:  filename (string), location (structure to be read in by the get
% sza functions.
%
% Output: time_info, line4, line4, sp (all matrices with rows corresponding
% with individual spectra and columsn corresponding with entries.

    % init the vectors
    time_info_tmp = [];
    spec_info = [];
    sp = [];
    is_sat = [];

    %opens the file and assigning a value to fid and check to see if a file
    %is found.  If not, continue to the next hour.
    fid = fopen(filename,'r');
    if fid == -1 %check to see if the file is found
        %disp(['No file: ' filename])
        return
    end
    
    % display filename
    disp(['Reading: ' filename])
    
    %skips first two headers
    fgetl(fid); fgetl(fid);
    
    line_nbr = 2;

    
    nbr_catches = 0;

    % now read through the file until we reach the end of the file
    while ~feof(fid)

        % read in header, remove ", and split it up into substrings by comma
        header_str = fgetl(fid);
        header_str = regexprep(header_str, '"', '');
        header = regexp(header_str, ',', 'split');
        
        if isempty(cell2mat(header)), continue, end
        
        if feof(fid), continue, end
        
        % read in the spectrum
        if version > 0 && version < 90
            sp_tmp = fscanf(fid, '%f', nbr_pix)';
            % and read in last empty line
            fgetl(fid);
            indy = find(sp_tmp > sat_lev);
            if length(indy) > 3;
                is_sat = [is_sat; true];
            else
                is_sat = [is_sat; false];
            end
            delta_line_nbr = nbr_pix+1;
        elseif version == 99
            sp_tmp_str = [];
            for i = 1:10
                sp_tmp_str = [sp_tmp_str ',' fgetl(fid)];
            end
            cell_tmp = regexp(sp_tmp_str,',','split');
            for w = 2:length(cell_tmp)
                sp_tmp(w) = str2num(cell_tmp{w});
            end
            nbr_pix = length(sp_tmp);
            indy = find(sp_tmp > sat_lev);
            if length(indy) > 3;
                is_sat = [is_sat; true];
            else
                is_sat = [is_sat; false];
            end
            delta_line_nbr = 11;
        elseif version == 0
            file_format = '%f';
            for i = 2:nbr_cols
                file_format = [file_format ',%f'];
            end
            sp_tmp=fscanf(fid, file_format, [nbr_cols,nbr_pix]);    
            [indx, indy] = find(sp_tmp > sat_lev);
            indx = unique(indx);
            sp_tmp(indx,:) = 0;
            sp_tmp = sum(sp_tmp);
            if length(indy > 3)
                is_sat = [is_sat; true];
            else
                is_sat = [is_sat; false];
            end
            delta_line_nbr = nbr_pix+1;
        end
        
        if isempty(sp_tmp); continue; end
        
        format_flag = 0;
        if length(sp_tmp) == 1
            format_flag = 1;
        elseif length(sp_tmp) < nbr_pix
            format_flag = 2;
        else
            sp_tmp = [sp_tmp, zeros(1,2048 - nbr_pix)];
            try
                st_date = cell2mat(header(1));
                if length(header) < 17
                    end_date = st_date;
                else
                    end_date = cell2mat(header(3));
                end
                if length(st_date) < 10
                    if version == 99
                        st_date = [st_date(1:6) '19' st_date(7:8)];
                        end_date = [end_date(1:6) '19' end_date(7:8)];
                    else
                        st_date =[st_date(1:6) '20' st_date(7:8)];
                        end_date =[end_date(1:6) '20' end_date(7:8)];
                    end
                    st_date_vec = [str2num(st_date(7:10)) ...
                        str2num(st_date(4:5)) str2num(st_date(1:2)) 0 0 0]; 
                    end_date_vec = [str2num(end_date(7:10)) ...
                        str2num(end_date(4:5)) str2num(end_date(1:2)) 0 0 0]; 
                else
                    st_date_vec = datevec(st_date, date_format);
                    end_date_vec = datevec(end_date, date_format);
                end
          
                sp = [sp; sp_tmp];
            catch
                format_flag = 1;
            end
        end
        
        flag_v0 = 0;
        if format_flag > 0
            nbr_catches = nbr_catches + 1;
            
            if nbr_catches > 4,
                disp('PROBLEM WITH FILE!')
                fclose(fid);
                return
            end
            
            disp('Encountered different file format')
            
            if nbr_catches == 2
                disp('Resetting number of spectral columns')
                fgetl(fid);
                sp_line = fgetl(fid);
                ind = findstr(',',sp_line);
                nbr_cols = length(ind)+1;
                flag_v0 = 1;
            end
            
            fclose(fid);
            fid = fopen(filename,'r');
            for i = 1:line_nbr; fgetl(fid); end
            if version == 0,
                if flag_v0 == 0
                    version = 1;
                end
                continue
            elseif version == 1,
                version = 0;
                continue
            else
                break
            end
        end
        
        nbr_catches = 0;
       
        tmp = datevec(header(2));%separate time
        st_date_vec(4:6) = tmp(4:6);
        
        if length(header) < 17
            end_date_vec(4:6) = st_date_vec(4:6);
        else
            tmp = datevec(header(4));%separate time
            end_date_vec(4:6) = tmp(4:6);
        end

        % now make line vectors corresponding to future L1 header lines
        % note that the date is flipped compared to the standard date vec
        % ie: dd, mm, yyyy, hh, mn, ss
        time_info_tmp = [time_info_tmp;...
            st_date_vec(1:3), st_date_vec(4:6), ...
            end_date_vec(1:3), end_date_vec(4:6)];
        
        % read in the fifth line of the L1 header
        tmp = [];
        if version == 99 && length(header) < 17
            for i = 3:15, tmp = [tmp str2num(cell2mat(header(i)))];  end
        else
            for i = 5:17, tmp = [tmp str2num(cell2mat(header(i)))];  end
        end
        spec_info = [spec_info; tmp]; % first index is shutter

        line_nbr = line_nbr + delta_line_nbr;
    end
    
    fclose(fid);
    
    time_info=table;
    time_info.start_time=datetime(time_info_tmp(:,1:6));
    time_info.end_time=datetime(time_info_tmp(:,7:12));
    
    spec_info=array2table(spec_info(:,4:12),'VariableNames',{'groove','turret','blaze',...
              'centre','t_int','no_accum','mean_T','min_T','max_T'});
          
    spec_info(:,2:3)=[];
    
end
