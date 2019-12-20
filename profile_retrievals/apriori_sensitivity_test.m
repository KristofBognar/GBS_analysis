
comp_1ppt=0;

plot_col=1; % 1: plot columns, 0: plot ratios

load('/home/kristof/work/profile_retrievals/profile_results/tracegas_profiles_all.mat')
times_orig=times;
info_orig=info;

% partial column height (top layer still included in the partial column)
top=4; % 600m layer, includes lab
% top=3; % 400m layer

part_prof_orig=prof_nd*20000; % layers are 200m thick
ratio_orig=sum(part_prof_orig(1:top,:))'./sum(part_prof_orig(top+1:end,:))';

if comp_1ppt
    %%
    
    figure, hold on
    if plot_col  
        plot(times_orig,info_orig.col,'ko')
    else
        plot(times_orig,ratio_orig,'ko')
    end
    disp('All surf. conc. from 5 ppt to 1 ppt')
    
    load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_surf_ppt_5to1/',...
          'tracegas/profiles_1_filt.mat'])
      
    if plot_col  
        plot(times,info.col,'rx')

        [~,inda,indb]=intersect(times,times_orig);

        abs=info.col(inda)-info_orig.col(indb);
        rel=((info.col(inda)./info_orig.col(indb))-1)*100;

        disp(mean(abs))
        disp(mean(rel))
    else
        part_prof=prof_nd*20000; 
        ratio=sum(part_prof(1:top,:))'./sum(part_prof(top+1:end,:))';

        plot(times,ratio,'rx')
    end

else
    %% scale height

    figure, hold on

    if plot_col  
        plot(times_orig,info_orig.col,'ko')
    else
        plot(times_orig,ratio_orig,'ko')
    end

    disp('Scale height')


    %%
    load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_scale_height_0p5/',...
          'tracegas/profiles_1_filt.mat'])
    if plot_col  
        plot(times,info.col,'rx')

        [~,inda,indb]=intersect(times,times_orig);

        abs=info.col(inda)-info_orig.col(indb);
        rel=((info.col(inda)./info_orig.col(indb))-1)*100;

        disp(mean(abs))
        disp(mean(rel))
    else
        part_prof=prof_nd*20000; 
        ratio=sum(part_prof(1:top,:))'./sum(part_prof(top+1:end,:))';

        plot(times,ratio,'rx')
    end

    %%
    load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_scale_height_1/',...
          'tracegas/profiles_1_filt.mat'])
    if plot_col  
        plot(times,info.col,'b+')

        [~,inda,indb]=intersect(times,times_orig);

        abs=info.col(inda)-info_orig.col(indb);
        rel=((info.col(inda)./info_orig.col(indb))-1)*100;

        disp(mean(abs))
        disp(mean(rel))
    else
        part_prof=prof_nd*20000; 
        ratio=sum(part_prof(1:top,:))'./sum(part_prof(top+1:end,:))';

        plot(times,ratio,'b+')
    end

    %%
    load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_scale_height_1p5/',...
          'tracegas/profiles_1_filt.mat'])
    if plot_col  
        plot(times,info.col,'gx')

        [~,inda,indb]=intersect(times,times_orig);

        abs=info.col(inda)-info_orig.col(indb);
        rel=((info.col(inda)./info_orig.col(indb))-1)*100;

        disp(mean(abs))
        disp(mean(rel))
    else
        part_prof=prof_nd*20000; 
        ratio=sum(part_prof(1:top,:))'./sum(part_prof(top+1:end,:))';

        plot(times,ratio,'gx')
    end

    %%
    load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_scale_height_2/',...
          'tracegas/profiles_1_filt.mat'])
    if plot_col  
        plot(times,info.col,'c+')

        [~,inda,indb]=intersect(times,times_orig);

        abs=info.col(inda)-info_orig.col(indb);
        rel=((info.col(inda)./info_orig.col(indb))-1)*100;

        disp(mean(abs))
        disp(mean(rel))
    else
        part_prof=prof_nd*20000; 
        ratio=sum(part_prof(1:top,:))'./sum(part_prof(top+1:end,:))';

        plot(times,ratio,'c+')
    end


    legend('orig','h=0.5','h=1','h=1.5','h=2')

    set(findobj(gca,'type','line'),'linew',2)


    %% surface conc.

    figure, hold on
    if plot_col  
        plot(times_orig,info_orig.col,'ko')
    else
        plot(times_orig,ratio_orig,'ko')
    end

    disp('Surface concentration')

    %%
    load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_surf_ppt_1/',...
          'tracegas/profiles_1_filt.mat'])
    if plot_col  
        plot(times,info.col,'rx')

        [~,inda,indb]=intersect(times,times_orig);

        abs=info.col(inda)-info_orig.col(indb);
        rel=((info.col(inda)./info_orig.col(indb))-1)*100;

        disp(mean(abs))
        disp(mean(rel))
    else
        part_prof=prof_nd*20000; 
        ratio=sum(part_prof(1:top,:))'./sum(part_prof(top+1:end,:))';

        plot(times,ratio,'rx')
    end

    %%
    load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_surf_ppt_5/',...
          'tracegas/profiles_1_filt.mat'])
    if plot_col  
        plot(times,info.col,'b+')

        [~,inda,indb]=intersect(times,times_orig);

        abs=info.col(inda)-info_orig.col(indb);
        rel=((info.col(inda)./info_orig.col(indb))-1)*100;

        disp(mean(abs))
        disp(mean(rel))
    else
        part_prof=prof_nd*20000; 
        ratio=sum(part_prof(1:top,:))'./sum(part_prof(top+1:end,:))';

        plot(times,ratio,'b+')
    end

    %%
    load(['/home/kristof/work/profile_retrievals/profile_results/eureka_1_surf_ppt_10/',...
          'tracegas/profiles_1_filt.mat'])
    if plot_col  
        plot(times,info.col,'gx')

        [~,inda,indb]=intersect(times,times_orig);

        abs=info.col(inda)-info_orig.col(indb);
        rel=((info.col(inda)./info_orig.col(indb))-1)*100;

        disp(mean(abs))
        disp(mean(rel))
    else
        part_prof=prof_nd*20000; 
        ratio=sum(part_prof(1:top,:))'./sum(part_prof(top+1:end,:))';

        plot(times,ratio,'gx')
    end



    legend('orig','sc=1','sc=5','sc=10')

    set(findobj(gca,'type','line'),'linew',2)

end
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

