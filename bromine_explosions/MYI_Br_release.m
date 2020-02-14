% various tests for justifying hypotheses of Br release from multi-year ice

load('/home/kristof/work/BEEs/BEE_dataset_all.mat')
bee_dataset(bee_dataset.times.Year==2015,:)=[];


frac_myi=bee_dataset.MYSI_5day./(bee_dataset.FYSI_5day+bee_dataset.MYSI_5day);

% ind_myi=(frac_myi>0.9 & bee_dataset.N_SE_rest==1);
ind_myi=(frac_myi>0.9);
ind_fyi=(frac_myi<0.1);

figure

ind_tmp=(ind_myi & bee_dataset.o3_surf<=15);
sizes=bee_dataset.MYSI_5day(ind_tmp)/max(bee_dataset.MYSI_5day(ind_tmp));
sizes=(sizes*100)+5;

scatter3(bee_dataset.mixing_height_5day(ind_tmp),...
         bee_dataset.frac_in_mix_5day(ind_tmp),...
         bee_dataset.bro_col(ind_tmp),sizes,'r'); hold on
  

ind_tmp=(ind_myi & bee_dataset.o3_surf>15);
sizes=bee_dataset.MYSI_5day(ind_tmp)/max(bee_dataset.MYSI_5day(ind_tmp));
sizes=(sizes*100)+5;

scatter3(bee_dataset.mixing_height_5day(ind_tmp),...
         bee_dataset.frac_in_mix_5day(ind_tmp),...
         bee_dataset.bro_col(ind_tmp),sizes,'k');
  
grid on
xlabel('hmix')
ylabel('fmix')
zlabel('BrO')











