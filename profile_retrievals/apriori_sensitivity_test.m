
load('/home/kristof/work/profile_retrievals/profile_results/tracegas_profiles_all.mat')
times_orig=times;
info_orig=info;

%% scale height

figure, hold on
plot(times_orig,info_orig.col,'ko')

disp('Scale height')

%%
load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_scale_height_0p5/',...
      'tracegas/profiles_1_filt.mat'])
plot(times,info.col,'rx')

[~,inda,indb]=intersect(times,times_orig);

abs=info.col(inda)-info_orig.col(indb);
rel=((info.col(inda)./info_orig.col(indb))-1)*100;

disp(mean(abs))
disp(mean(rel))

%%
load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_scale_height_1/',...
      'tracegas/profiles_1_filt.mat'])
plot(times,info.col,'b+')

[~,inda,indb]=intersect(times,times_orig);

abs=info.col(inda)-info_orig.col(indb);
rel=((info.col(inda)./info_orig.col(indb))-1)*100;

disp(mean(abs))
disp(mean(rel))

%%
load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_scale_height_1p5/',...
      'tracegas/profiles_1_filt.mat'])
plot(times,info.col,'gx')

[~,inda,indb]=intersect(times,times_orig);

abs=info.col(inda)-info_orig.col(indb);
rel=((info.col(inda)./info_orig.col(indb))-1)*100;

disp(mean(abs))
disp(mean(rel))

%%
load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_scale_height_2/',...
      'tracegas/profiles_1_filt.mat'])
plot(times,info.col,'c+')

[~,inda,indb]=intersect(times,times_orig);

abs=info.col(inda)-info_orig.col(indb);
rel=((info.col(inda)./info_orig.col(indb))-1)*100;

disp(mean(abs))
disp(mean(rel))


legend('orig','h=0.5','h=1','h=1.5','h=2')

set(findobj(gca,'type','line'),'linew',2)


%% surface conc.

figure, hold on
plot(times_orig,info_orig.col,'ko')

disp('Surface concentration')

%%
load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_surf_ppt_1/',...
      'tracegas/profiles_1_filt.mat'])
plot(times,info.col,'rx')

[~,inda,indb]=intersect(times,times_orig);

abs=info.col(inda)-info_orig.col(indb);
rel=((info.col(inda)./info_orig.col(indb))-1)*100;

disp(mean(abs))
disp(mean(rel))

%%
load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_surf_ppt_5/',...
      'tracegas/profiles_1_filt.mat'])
plot(times,info.col,'b+')

[~,inda,indb]=intersect(times,times_orig);

abs=info.col(inda)-info_orig.col(indb);
rel=((info.col(inda)./info_orig.col(indb))-1)*100;

disp(mean(abs))
disp(mean(rel))

%%
load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_surf_ppt_10/',...
      'tracegas/profiles_1_filt.mat'])
plot(times,info.col,'gx')

[~,inda,indb]=intersect(times,times_orig);

abs=info.col(inda)-info_orig.col(indb);
rel=((info.col(inda)./info_orig.col(indb))-1)*100;

disp(mean(abs))
disp(mean(rel))



legend('orig','sc=1','sc=5','sc=10')

set(findobj(gca,'type','line'),'linew',2)

%% find good test days
% % load('/home/kristof/work/documents/paper_bro/data/BEE_dataset_all.mat')
% % bee_dataset(bee_dataset.times.Year==2015,:)=[];
% % apriori=get_BrO_apriori(bee_dataset.times);
% % 
% % dt_ind=apriori.scale_h;
% % dt_ind(dt_ind==0.5)=0;
% % 
% % surf_ind=apriori.surf_ppt;
% % surf_ind(surf_ind==1)=0;
% % surf_ind(surf_ind==5)=1;
% % 
% % % ind=(bee_dataset.sonde_dT>12);
% % % gscatter(bee_dataset.times(ind),bee_dataset.bro_col(ind),dt_ind(ind),'rb','oo')
% % 
% % ind=(1:length(bee_dataset.sonde_dT));
% % gscatter(bee_dataset.times(ind),bee_dataset.bro_col(ind),surf_ind(ind),'rb','oo')

