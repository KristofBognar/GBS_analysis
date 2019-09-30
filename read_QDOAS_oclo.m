function [dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS_oclo(file_nm, col, filt,version,source,save_fig,working_dir,real_sza)
% [dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS(file_nm, col, filt,2)
% version = 1 for DOAS data
% version = 2 for MAX-DOAS data, this function will be used as subfunction
% for read_QDOAS_maxdoas.m


% This function creates the daily data structure that the rest of the
% software uses to analyze VCDs.  

% for DOAS reading
% ex:[dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS('UE_2011_0660-70.asc',col_no2_p0, filt_good,1,0);
% ex: Ozone : [dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS('UE_2011_600gr_fx.dat', col_o3_2, filt_good,1,0);

% INPUT:
%   file_nm:  string file-name of QDOAS output
%   col: structure containing column information listed below
    %   col.tot_nbr:  total number of columns in QDOAS output file
    %   col.fd: column for day of year
    %   col.sza: column for SZA
    %   col.saa: column for solar azimuth
    %   col.elev: column containing viewing elevation
    %   col.dscd: column for DSCD
    %   col.err: column for DSCD error
    %   col.rms: column for RMS
    %   col.x : column for the X DSCD if applicable (used for ozone error
    %   characterization!).
    %   col.shift : column containing the wavelength shift information
    %   col.stretch : column containing the wavelength shift information
    %   col.ref_sza: column containing the reference sza information
    %   col.year: column containing measurement year
%   filt: structure containing filtering information (if field is empty, do
%   not filter by this value)
    %   filt.elev: pointing elevation to be filtered by
    %   filt.rms_vec: vector containing two columns corresponding with SZA and RMS
    %   value.  This represents a threshold for RMS at the various SZA.
    %   filt.max_sza: maximum allowed SZA in the analysis (often set to 92 for
    %            AMF range)
    %   filt.shift: maximum allowed shift between reference and twilight
    %   spectrum.
%
% OUTPUT:
%    qdoas_filt: the full QDOAS data matrix corresponding to the
%               the filtered QDOAS data.
%    qdoas_raw: raw QDOAS data matrix
%    dscd_S: this is an object that contains the fields
%               which are relevant to the Langley analysis.
%               It contains columns the following columns all of which
%               should have the same length!
%                   dscd_S.day: day in local time
%                   dscd_S.fd :fractional day (UTC)
%                   dscd_S.ampm: (am = 0, pm = 1)
%                   dscd_S.sza: SZA
%                   dscd_S.saa: solar azimuth
%                   dscd_S.rms : rms
%                   dscd_S.mol_dscd : focus molecule dscd
%                   dscd_S.err : focus molecule dscd fitting error
%                   dscd_S.x : fit of X cross-section
%                   dscd_S.shift: shift in spectrum versus reference (nm)
%                   dscd_S.stretch: strech in spectrum versus reference
%                   (nm)
%                   dscd_S.ref_sza: SZA of reference spectrum
%                   dscd_S.year: year of measurements


    if nargin == 5 % the default setting won't save the figure
        save_fig = 0; % save the figures, 1 = yes, 0 = no
        working_dir = pwd; % this is the place you want save the figure
    end

    % Read in the QDOAS data
    fid = fopen(file_nm, 'r');
    fgetl(fid);
    
    if source=='C'
        % to read Cristen's files (2011 and before)
        qdoas_raw = (fscanf(fid,'%f', [45,inf]))';
        col.rms=13;
        col.ref_sza=12;
    else
        qdoas_raw = (fscanf(fid,'%f', [28,inf]))';
    end
    
    fclose(fid);
    
    % Now sort the data up by fractional day and ID whether there are
    % doubles of some values
    qdoas_raw = sortrows(qdoas_raw, col.fd);
    all_ind = 1:length(qdoas_raw(:,1));
    [a, unique_ind,b] = unique(qdoas_raw(:, col.fd));
    diff_ind = setdiff(all_ind, unique_ind);
    if ~isempty(diff_ind)
        disp('[WARNING]: File contains multiple entries taken at the same time.')
        disp('frac_day')
        disp('---------')
        disp(qdoas_raw(diff_ind, col.fd))
    end
    
    % filter out rms = 0, and dscd = 9999
    ind = find(qdoas_raw(:, col.rms) ~= 0 & ...
        qdoas_raw(:, col.oclo_dscd) ~= 9999 & ...
        qdoas_raw(:, col.oclo_err) ~= 9999 & ...
        qdoas_raw(:, col.sza) < filt.sza_max);
    qdoas_tmp = qdoas_raw(ind,:);
    
    % and filter by viewing elevation
    try
        ind = find( qdoas_tmp(:, col.elev) == filt.elev |...
            qdoas_tmp(:, col.elev) < -90);
    catch
        disp('No viwing elevation filter.')
        ind = 1:length(qdoas_tmp(:,1));
    end    
    qdoas_tmp = qdoas_tmp(ind, :);
    
    len_init = length(ind);
    
    % This block of the code deals with filtering the qdoas dataset
    % Filter by RMS?
    if version==1
        try
            ind = [];
            prev_sza = 0;
            for i = 1:length(filt.rms_vec(:,1))
                if i == 1,
                    ind = [ind; find(qdoas_tmp(:, col.sza) <= filt.rms_vec(i,1)...
                        & qdoas_tmp(:, col.sza) >= prev_sza...
                        & qdoas_tmp(:, col.rms) <= filt.rms_vec(i,2))];        
                        prev =filt.rms_vec(i,1);
                else
                    ind = [ind; find(qdoas_tmp(:, col.sza) <= filt.rms_vec(i,1)...
                        & qdoas_tmp(:, col.sza) > filt.rms_vec(i-1,1)...
                        & qdoas_tmp(:, col.sza) >= prev_sza...
                        & qdoas_tmp(:, col.sza) < filt.sza_max...
                        & qdoas_tmp(:, col.rms) <= filt.rms_vec(i,2))];        
                        prev =filt.rms_vec(i,1);
                end

            end
            ind = sort(ind');
        catch
            disp('No RMS filter.  Either filt.rms_vec or filt.sza_max are empty!')
            ind = 1:length(qdoas_tmp(:,1));
        end
    elseif version==2
        ind = find(qdoas_tmp(:, col.rms) <= filt.rms);
    else
        ind = 1:length(qdoas_tmp(:,1));
    end
    temp = qdoas_tmp(ind, :);
    
    % and filter by shift
    try
        ind = find( abs(temp(:, col.shift)) < filt.shift);
        temp = temp(ind, :);
    catch
        disp('No shift filter.')
        ind = 1:length(temp(:,1));
    end
    qdoas_filt = temp;
    
    if version ==1
        % See if user passed actual SZA filter, plot accordingly
        if nargin==7
            sza_ind=find(qdoas_filt(:,col.sza) <= real_sza(2) ...
                & qdoas_filt(:,col.sza) >= real_sza(1));
        else
            sza_ind=1:size(qdoas_filt,1);
        end    

        % Make figures showing the results of the filtering
        figure
        subplot(3,1,1)
        hold on
        plot(qdoas_raw(:,col.fd), qdoas_raw(:, col.oclo_dscd), '.', 'color','r')
        plot(qdoas_filt(sza_ind,col.fd), qdoas_filt(sza_ind, col.oclo_dscd), 'b.')
        legend('Raw', 'Filt')
        xlabel('day')
        ylabel('DSCD (mol/cm^2)')
        subplot(3,1,2)
        hold on
        plot(qdoas_raw(:,col.fd), qdoas_raw(:, col.rms), '.', 'color','r')
        plot(qdoas_filt(sza_ind,col.fd), qdoas_filt(sza_ind, col.rms), 'b.')
        legend('Raw', 'Filt')
        xlabel('day')
        ylabel('RMS')
        ylim([0,0.005])
        subplot(3,1,3)
        hold on
        plot(qdoas_raw(:,col.fd), qdoas_raw(:, col.sza), '.', 'color','r')
        plot(qdoas_filt(sza_ind,col.fd), qdoas_filt(sza_ind, col.sza), 'b.')
        legend('Raw', 'Filt')
        xlabel('day')
        ylabel('SZA')
        cd(working_dir);
        print_setting(1/2,save_fig,['results_of_filtering']);
    end
    % now sort data by day
    
    % sort data by azimuth and fractional day
    qdoas_filt = sortrows(qdoas_filt, col.saa);
    qdoas_filt = sortrows(qdoas_filt, col.fd);
    % find indices of data that should be on previous day
    prev_day_i = find(qdoas_filt(:, col.saa) > 0 ...
        & rem(qdoas_filt(:, col.fd),1) < 0.25);
    
    % find indices of data that should be pm
    pm_i = find(qdoas_filt(:, col.saa) > 0);
    
    % now create new data columns: 
    % local time day, fractional day, am/pm, sza, az, rms, dscd, err
    L = length(qdoas_filt(:,1));
    dscd_S.day = floor(qdoas_filt(:, col.fd));
    dscd_S.fd = qdoas_filt(:,col.fd); 
    dscd_S.ampm = zeros(L, 1);
    dscd_S.sza = qdoas_filt(:, col.sza);
    dscd_S.saa = qdoas_filt(:, col.saa);
    dscd_S.rms = qdoas_filt(:,col.rms);
    dscd_S.mol_dscd = qdoas_filt(:, col.oclo_dscd);
    dscd_S.err = qdoas_filt(:, col.oclo_err);
    if version == 2 %for MAX-DOAS O4 reading
        dscd_S.mol_dscd_o4 = qdoas_filt(:, col.dscd_o4);
        dscd_S.err_o4 = qdoas_filt(:,col.err_o4);
%         % for MAX-DOAS O3 reading
%         dscd_S.mol_dscd_o3 = qdoas_filt(:, col.dscd_o3);
%         dscd_S.err_o3 = qdoas_filt(:,col.err_o3);
%         % for MAX-DOAS NO2 reading
%         dscd_S.mol_dscd_no2 = qdoas_filt(:, col.dscd_no2);
%         dscd_S.err_no2 = qdoas_filt(:,col.err_no2);
%         % for MAX-DOAS OClO reading
%         dscd_S.mol_dscd_oclo = qdoas_filt(:, col.dscd_oclo);
%         dscd_S.err_oclo = qdoas_filt(:,col.err_oclo);
%         % for MAX-DOAS Ring reading
%         dscd_S.mol_dscd_ring = qdoas_filt(:, col.dscd_ring);
%         dscd_S.err_ring = qdoas_filt(:,col.err_ring);
    end;
    
    try
        dscd_S.x = qdoas_filt(:, col.x);
    catch
        disp('No column for x cross-section supplied')
    end
    try
        dscd_S.shift = qdoas_filt(:, col.shift);
    catch
        disp('No column for shift supplied')
    end
    try
        dscd_S.stretch = qdoas_filt(:, col.stretch);
    catch
        disp('No column for stretch supplied')
    end
    try
        dscd_S.ref_sza = qdoas_filt(:, col.ref_sza);
    catch
        disp('No column for ref SZA supplied')
    end
    try
        dscd_S.year = qdoas_filt(:, col.year);
    catch
        disp('No column for year supplied')
    end
    dscd_S.ampm(pm_i) = 1;
    dscd_S.day(prev_day_i) = dscd_S.day(prev_day_i) - 1;
    
    if version ==1
        % make figure
        figure
        subplot(2,1,1)
        hold on
        i1 = find(dscd_S.ampm == 0);
        plot(dscd_S.fd(i1), dscd_S.mol_dscd(i1),'.','color','b') 
        i2 = find(dscd_S.ampm == 1);
        plot(dscd_S.fd(i2), dscd_S.mol_dscd(i2),'.','color','r') 
        xlabel('Fractional day')
        ylabel('DSCD (mol/cm^2)')
        legend('am','pm')
        subplot(2,1,2)
        hold on
        plot(dscd_S.fd(i1), dscd_S.sza(i1),'.','color','b') 
        plot(dscd_S.fd(i2), dscd_S.sza(i2),'.','color','r')
        xlabel('Fractional day')
        ylabel('SZA')
        legend('am','pm')
        print_setting(1/2,save_fig,['filtered_DSCDs']);
        len_fin = length(dscd_S.sza);

        disp(['Initial # points (after elev and SZA max filter): ' num2str(len_init)])
        disp(['Final # points (after all filters): '  num2str(len_fin)])
        disp([num2str( ((len_init - len_fin) / len_init) * 100, 5) '% of data filtered']);
    end
end