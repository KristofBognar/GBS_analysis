function data=merge_RD_QDOAS_output( instrument, year )
%merge_RD_QDOAS_output( year ) Merge QDOAS output from RD processing into
%single yearly file, to be used in yearly VCD retrieval
%
% instrument = 'UT-GBS' or 'PEARL-GBS'
% year = integer year number
% 
% input: RD files saved as .mat in folder specified in rd_dir
%
% output: yearly files saved as .mat in folder specified in yearly_dir
%

%% setup

if ismac
    rd_dir='/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/QDOAS_results/NDACC_RD_tables';
    yearly_dir='/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/QDOAS_results/yearly_tables';
elseif isunix
    rd_dir='/home/kristof/work/GBS/QDOAS_results/NDACC_RD_tables';
    yearly_dir='/home/kristof/work/GBS/QDOAS_results/yearly_tables';
end

data_all=[];

if strcmp(instrument,'PEARL-GBS')
    uv_str='_UV';
    disp('Assuming PEARL-GBS data is UV')
elseif strcmp(instrument,'UT-GBS')
    uv_str='';
    disp('Assuming UT-GBS data is VIS')
else
    error('Give correct instrument name!');
end

%% read files

cd(rd_dir)

tmp = dir([instrument '_' num2str(year) '*.mat']); 
f_list = {tmp.name}; % cell array of file names

for i=f_list
    
    load(i{1})
    data_all=[data_all; data];
    
end

data=sortrows(data_all,3);

save([yearly_dir '/' instrument '_' num2str(year) uv_str '.mat'],'data');

end

