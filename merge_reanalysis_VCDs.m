function merge_reanalysis_VCDs(instr, tg, year, is_rd, do_filter)
% merge_reanalysis_VCDs(instr, tg, year, is_rd, do_filter)
% merge and filter yearly VCD values

if nargin==4, do_filter=1; end

%% setup

%instrument
% instr='UT-GBS';
% instr='PEARL-GBS';

% trace gas (1: ozone, 2: NO2, 3: NO2 UV)
% tg=2;
tgstr={'O3','NO2','NO2_UV'};

% year
% year='2019'; % set to empty string to include all files

% filters
% do_filter=true;

addRCDerr=true;
do_cams_filter=true;

%% output

% VCD directory
if ismac
    if ~is_rd
        vcd_dir='/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/VCD_results/';
    else
        vcd_dir='/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/VCD_results/NDACC_RD/';
    end
elseif isunix
    if ~is_rd
        vcd_dir='/home/kristof/work/GBS/VCD_results/';
    else
        vcd_dir='/home/kristof/work/GBS/VCD_results/NDACC_RD/';
    end
end
    

% output file name (remove file if already exists)

if do_filter % regular merged files, with bad values filtered out
    fname=[vcd_dir instr '_' tgstr{tg} '_VCD_' year 'all.mat'];
else % unfiltered data for DMP input
    fname=[vcd_dir instr '_' tgstr{tg} '_VCD_' year '_unfiltered.mat'];
end
if exist(fname,'file'), delete(fname), end

%% find files
reanalysis=[];

% make list of VCD files
tmp = dir([vcd_dir instr '_' tgstr{tg} '_VCD_' year '*.mat']); 
f_list = {tmp.name}; % cell array of file names

if strcmp(instr,'UT-GBS')
    instr_in=1;
else
    instr_in=2;
end

%% loop over VCD files
cd(vcd_dir)

for file=f_list
    
    % skip merged files if they are included in the list (of files to be merged)
    if ~isempty(strfind(file{1},'all')), continue, end
    if ~isempty(strfind(file{1},'unfiltered')), continue, end

    load(file{1})

    if do_filter
        % filter out NaNs and bad RCD values
        [ind_goodvcd,VCD_table2] = filter_VCD_output( tg, VCD_table, rcd_S, instr_in,...
                                                     addRCDerr, do_cams_filter);
    
        VCD_filt=VCD_table2(ind_goodvcd,:);
        
    else
        
        VCD_filt=VCD_table;
    
    end
    
%     reanalysis=vertcat(reanalysis,VCD_filt);
    reanalysis=[reanalysis;VCD_filt];

end

%% add time information
reanalysis.mjd2k=ft_to_mjd2k(reanalysis.fd-1,reanalysis.year);
reanalysis.fractional_time=reanalysis.fd-1;

%% save file
reanalysis=sortrows(reanalysis,'mjd2k');
save(fname,'reanalysis');

clearvars
