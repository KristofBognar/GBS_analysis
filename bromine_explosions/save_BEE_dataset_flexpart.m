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
bro_prctile=[];
for i=10:10:90
    bro_prctile(i/10)=prctile(bee_dataset.bro_col,i);
end

bro_mean=mean(bee_dataset.bro_col);

% get dp>0.5 micron aer concentration percentiles
ssa_prctile=[];
for i=10:10:90
    ssa_prctile(i/10)=prctile(bee_dataset.aer_halfmicron,i);
end

ssa_mean=mean(bee_dataset.aer_halfmicron);

% get supermicron aer concentration percentiles
ssa_prctile_sm=[];
for i=10:10:90
    ssa_prctile_sm(i/10)=prctile(bee_dataset.aer_supermicron,i);
end

ssa_mean_sm=mean(bee_dataset.aer_supermicron);

%% group dataset

% group by wind direction (take direction that's there for over 50% of
% the profiles
[ind_wdir]=group_BrO(bee_dataset, run_start, run_end, 'wdir', 0.5);
ind_wdir(ind_wdir==0)=3; % merge other>50% and none>50%

% get mean wind speed
[~,mean_wspd]=group_BrO(bee_dataset, run_start, run_end, 'wspd',8);

%%% group by BrO column
% percentiles
[ind_bro_pc,run_mean_bro]=group_BrO(bee_dataset, run_start, run_end, 'bro', bro_prctile);
% mean
ind_bro_m=group_BrO(bee_dataset, run_start, run_end, 'bro', bro_mean);

%%% group by Dp>0.5 micron aer concentration
% percentiles
[ind_ssa_pc,run_mean_ssa]=group_BrO(bee_dataset, run_start, run_end, 'ssa_hm', ssa_prctile);
% mean
ind_ssa_m=group_BrO(bee_dataset, run_start, run_end, 'ssa_hm', ssa_mean);

%%% group by supermicron aer concentration
% percentiles
[ind_ssa_pc_sm,run_mean_ssa_sm]=group_BrO(bee_dataset,run_start,run_end,'ssa_sm',ssa_prctile_sm);
% mean
ind_ssa_m_sm=group_BrO(bee_dataset, run_start, run_end, 'ssa_sm', ssa_mean_sm);

%%% group by surface ozone
[ind_o3,run_mean_o3]=group_BrO(bee_dataset, run_start, run_end, 'o3', [5,10,20]);

%%% group by surface to lab T inversion strength
[ind_dT,run_mean_dT]=group_BrO(bee_dataset, run_start, run_end, 'inv', 7);


%% save data
 
bee_fp=table();

bee_fp.mean_time=run_times';
bee_fp.wdir=ind_wdir;
bee_fp.wspd=mean_wspd;
bee_fp.bro_pc=ind_bro_pc;
bee_fp.bro_m=ind_bro_m;
bee_fp.bro_mean_col=run_mean_bro;
bee_fp.ssa_pc=ind_ssa_pc;
bee_fp.ssa_m=ind_ssa_m;
bee_fp.ssa_mean=run_mean_ssa;
bee_fp.ssa_pc_sm=ind_ssa_pc_sm;
bee_fp.ssa_m_sm=ind_ssa_m_sm;
bee_fp.ssa_mean_sm=run_mean_ssa_sm;
bee_fp.o3=ind_o3;
bee_fp.o3_mean=run_mean_o3;
bee_fp.dT=ind_dT;
bee_fp.dT_mean=run_mean_dT;

%% add ice /water/land contact times

for f=1:5 % back traj length in days

    if f==3
        
        load(['/home/kristof/work/BEEs/flexpart_SI_contact/FP_FYSI_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.fysi_' num2str(f) 'day=FP_SI_contact.contact;']);

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/FP_MYSI_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.mysi_' num2str(f) 'day=FP_SI_contact.contact;']);

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/approximate/FP_water_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.approx_water_' num2str(f) 'day=FP_SI_contact.contact;']);

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/approximate/FP_land_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.approx_land_' num2str(f) 'day=FP_SI_contact.contact;']);
        
    else
        
        load(['/home/kristof/work/BEEs/flexpart_SI_contact/approximate/FP_FYSI_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.approx_fysi_' num2str(f) 'day=FP_SI_contact.contact;']);

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/approximate/FP_MYSI_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.approx_mysi_' num2str(f) 'day=FP_SI_contact.contact;']);

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/approximate/FP_water_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.approx_water_' num2str(f) 'day=FP_SI_contact.contact;']);

        load(['/home/kristof/work/BEEs/flexpart_SI_contact/approximate/FP_land_contact_'...
              num2str(f) 'day.mat'])
        eval(['bee_fp.approx_land_' num2str(f) 'day=FP_SI_contact.contact;']);
    
    end
end

%% save file

% remove NaNs (FP rauns where all BrO measurements were filtered out)
% bee_fp(isnan(bee_fp.bro_mean_col),:)=[];

save('/home/kristof/work/BEEs/BEE_dataset_flexpart.mat','bee_fp');

end

