function save_BEE_dataset_flexpart()
% save_BEE_dataset_flexpart(): Group BEE dataset by various parameters
% for time periods coresponding to FLEXPART back run start times
%
% Refer to this file and group_BrO.m to figure out what the hell the output
% means...


%% setup

flex_dir='/home/kristof/atmosp_servers/export/data/home/kbognar/FLEXPART_10.02/';
flex_folder='BrO_back_runs_v1';


% load BrO dataset
load('/home/kristof/work/BEEs/BEE_dataset_all.mat')

ymin=min(year(bee_dataset.times));
ymax=max(year(bee_dataset.times));

% cd([flex_dir flex_folder]);

% read BrO times corresponding to flecpart runs
try
    load([flex_dir flex_folder '/flexpart_times_' num2str(ymin) '-' num2str(ymax) '.mat']);
catch
    error('Merge/split flexpart times into one file')
end
    
%% get info to group dataset by

% get BrO column percentiles
bro_25=prctile(bee_dataset.bro_col,25);
bro_50=prctile(bee_dataset.bro_col,50);
bro_75=prctile(bee_dataset.bro_col,75);
bro_90=prctile(bee_dataset.bro_col,90);
bro_mean=mean(bee_dataset.bro_col);

% get supermicron aer concentration percentiles
ssa_25=prctile(bee_dataset.OPC_supermicron,25);
ssa_50=prctile(bee_dataset.OPC_supermicron,50);
ssa_75=prctile(bee_dataset.OPC_supermicron,75);
ssa_90=prctile(bee_dataset.OPC_supermicron,90);
ssa_mean=mean(bee_dataset.OPC_supermicron);

%% group dataset

% group by wind direction (take direction that's there for over 50% of
% the profiles
[ind_wdir]=group_BrO(bee_dataset, run_start, run_end, 'wdir', 0.5);
ind_wdir(ind_wdir==0)=3; % merge other>50% and none>50%

%%% group by BrO column
% percentiles
[ind_bro_pc,run_mean_bro]=group_BrO(bee_dataset, run_start, run_end, 'bro', [bro_25, bro_50, bro_75, bro_90]);
% mean
ind_bro_m=group_BrO(bee_dataset, run_start, run_end, 'bro', bro_mean);

%%% group by supermicron aer concentration
% percentiles
[ind_ssa_pc,run_mean_ssa]=group_BrO(bee_dataset, run_start, run_end, 'ssa', [ssa_25, ssa_50, ssa_75, ssa_90]);
% mean
ind_ssa_m=group_BrO(bee_dataset, run_start, run_end, 'ssa', ssa_mean);

%%% group by surface ozone
[ind_o3,run_mean_o3]=group_BrO(bee_dataset, run_start, run_end, 'o3', [5,10,20]);

%%% group by surface to lab T inversion strength
[ind_dT,run_mean_dT]=group_BrO(bee_dataset, run_start, run_end, 'inv', 7);

%% save data
 
bee_fp=table();

bee_fp.mean_time=run_times';
bee_fp.wdir=ind_wdir;
bee_fp.bro_pc=ind_bro_pc;
bee_fp.bro_m=ind_bro_m;
bee_fp.bro_mean_col=run_mean_bro;
bee_fp.ssa_pc=ind_ssa_pc;
bee_fp.ssa_m=ind_ssa_m;
bee_fp.ssa_mean=run_mean_ssa;
bee_fp.o3=ind_o3;
bee_fp.o3_mean=run_mean_o3;
bee_fp.dT=ind_dT;
bee_fp.dT_mean=run_mean_dT;


save('/home/kristof/work/BEEs/BEE_dataset_flexpart.mat','bee_fp');

end

