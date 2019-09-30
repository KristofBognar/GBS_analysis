
resid_files={'o3_res_101.ASC','o3_res_121.ASC','o3_res_141.ASC',...
             'o3_res_161.ASC','o3_res_181.ASC'};
         
days={'d101','d121','d141','d161','d181'};            

% average all residuals on only highest SZA
do_window=false;

% size of sza window (degrees)
window=5;

res_all=[];

for jj=1:length(resid_files)         
         
    % load residual file
    data=load(resid_files{jj});

    % assign variables
    sza=data(2:end,2);
    fd=data(2:end,3);

    % find not NaN entries (file contains full wavelength range outside fitting window)
    good_ind=find(~isnan(data(2,:)));
    good_ind(1:3)=[];

    lambda=data(1,good_ind);
    resid=data(2:end,good_ind);

    % average residuals
    if ~do_window
        mean_res=mean(resid,1);
        ind=1:length(fd);
    else
        ind=find(sza>=max(sza)-window);
        mean_res_5deg=mean(resid(ind,:),1);
    end

    if jj==1, lambda_all=lambda; res_all=zeros(size(lambda)); end

    figure(1)
    
    subplot(3,2,jj)

    for i=1:length(fd)
        plot(sza(i),mean((mean(resid(ind,:),1)./resid(i,:)) -1),'kx'), hold on
    end
    
    title(days{jj})
    xlabel('SZA (deg)')
    ylabel('\Delta residual (%)')
    ylim([-40,40])
    xlim([55,95])
 
    figure(2)
    
    subplot(3,2,jj)
    if ~do_window
        plot(lambda,mean_res), hold on
    else
        plot(lambda,mean_res_5deg), hold on
    end
    title(days{jj})
    xlabel('Lambda (nm)')
    ylabel('Residual')
    ylim([-0.01,0.011])
    
    out=[lambda',mean_res'];
    dlmwrite(['X_U1_2016_' days{jj} '_o3.xs'],out,'\t');
    
    res_all=res_all+interp1(lambda,mean_res,lambda_all);
        
end

res_all=res_all/length(resid_files);

ind=find(isnan(res_all));

lambda_all(ind)=[];
res_all(ind)=[];


out=[lambda_all',res_all'];
dlmwrite(['X_U1_2016_mean_o3.xs'],out,'\t');








