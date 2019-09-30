% merge yearly VCD values (don't filter or save)

error('use merge_reanalysis_VCDs.m')

%% control variables

%instrument
% instr='UT-GBS';
instr='PEARL-GBS';

year='2018';
% year='';

% trace gas (1: ozone, 2: NO2, 3: NO2 UV)
tg=3;
tgstr={'O3','NO2','NO2_UV'};

% VCD directory
% vcd_dir='/home/kristof/work/GBS/VCD_results/';
vcd_dir='/home/kristof/work/GBS/VCD_results/NDACC_RD/';

%% find files
reanalysis=[];

% make list of VCD files
temp = dir([vcd_dir instr '_' tgstr{tg} '_VCD_' year '*.mat']); 
f_list = {temp.name}; % cell array of file names

% remove filtered merged file
fname=[instr '_' tgstr{tg} '_VCD_all.mat'];
ind=find_in_cell(f_list,fname);
if ~isempty(ind), f_list(ind)=[]; end

%% loop over VCD files
cur_dir=pwd();
cd(vcd_dir)

for file=f_list
    
    load(file{1})
    
    reanalysis=[reanalysis;VCD_table];

end

%% add time information
reanalysis.mjd2k=ft_to_mjd2k(reanalysis.fd-1,reanalysis.year);
reanalysis.fractional_time=reanalysis.fd-1;

%% save file?

cd(cur_dir)
clearvars -except reanalysis
