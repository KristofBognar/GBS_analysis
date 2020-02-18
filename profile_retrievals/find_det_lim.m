% find approximate detection limit for BrO, using 10 days of results from
% August 2018 -- should contain minimal to no BrO

% load results
data_path='/home/kristof/work/profile_retrievals/profile_results/eureka_1_det_lim_test/';

% aerosol data
load([data_path 'aerosol/profiles_1_filt.mat'])
info_aer=info;
times_aer=times;

% tracegas data
load([data_path 'tracegas/profiles_1_filt.mat'])

% remove times when BrO retrieval failed
[times,ia,ib]=intersect(times,times_aer);

info=info(ia,:);
info_aer=info_aer(ib,:);

% filters (double check that it is the same as in save_BEE_dataset.m)
% exclude two highest datapoints: might actually be decently high BrO
ind_bad=find( info.DOFS<0.7 | info_aer.DOFS<0.7 | info_aer.col>5 | info.col>12*1e12);

times(ind_bad)=[];
info(ind_bad,:)=[];
info_aer(ind_bad,:)=[];

% fit gaussian curve to data
[y,edges] = histcounts(info.col,50);
x=(edges(1:end-1)+edges(2:end))/2;

[fitresult, gof] = fit( x', y', 'gauss1' );

det_lim=fitresult.b1+3*fitresult.c1;

figure
scatter(x,y), hold on
plot(fitresult)
plot([det_lim,det_lim],[0,max(y)])


disp(['Approx. BrO detection limit:'])
disp(det_lim)



