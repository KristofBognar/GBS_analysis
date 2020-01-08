% merge and filter yearly VCD values

%% control variables

%instrument
% instr='UT-GBS';
instr='PEARL-GBS';

% trace gas (1: ozone, 2: NO2, 3: NO2 UV)
tg=3;
tgstr={'O3','NO2','NO2_UV'};

% year
% year='2019'; % set to empty string to include all files
year=''; % set to empty string to include all files

% filters
do_filter=true;

addRCDerr=true;
do_cams_filter=true;

% VCD directory
% vcd_dir='/home/kristof/work/GBS/VCD_results/UT-GBS_reanalysis_old_err_budget/';
% vcd_dir='/home/kristof/work/GBS/VCD_results/PEARL-GBS_reanalysis_old_err_budget/';
vcd_dir='/home/kristof/work/GBS/VCD_results/';
% vcd_dir='/home/kristof/work/GBS/VCD_results/NDACC_RD/';

% output file name (remove file if already exists)
fname=[vcd_dir instr '_' tgstr{tg} '_VCD_' year 'all.mat'];
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

save(fname,'reanalysis');

clearvars
