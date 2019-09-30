function merge_GBSs()
%merge_GBSs merge VIS o3, no2, and UV no2 datasets from GBSs
%   

merge_o3=0;
merge_no2=0;
merge_no2uv=1;

% directories
vcd_dir='/home/kristof/work/GBS/VCD_results/';
save_dir='/home/kristof/work/satellite_validation/';

%% merge VIS ozone

if merge_o3
    try
        load([vcd_dir 'UT-GBS_O3_VCD_all.mat'])
        ut=reanalysis;

        load([vcd_dir 'PEARL-GBS_O3_VCD_all.mat'])
        p=reanalysis;

    catch
        error('One (or both) VIS O3 files missing')
    end

    gbs_o3=merge(ut,p);

    save([save_dir 'GBS_O3.mat'], 'gbs_o3')
end

%% merge VIS no2

if merge_no2
    try
        load([vcd_dir 'UT-GBS_NO2_VCD_all.mat'])
        ut=reanalysis;

        load([vcd_dir 'PEARL-GBS_NO2_VCD_all.mat'])
        p=reanalysis;

    catch
        error('One (or both) VIS NO2 files missing')
    end

    gbs_no2=merge(ut,p);

    save([save_dir 'GBS_NO2.mat'], 'gbs_no2')
end

%% merge UV no2

if merge_no2uv
    try
        load([vcd_dir 'UT-GBS_NO2_UV_VCD_all.mat'])
        ut=reanalysis;

        load([vcd_dir 'PEARL-GBS_NO2_UV_VCD_all.mat'])
        p=reanalysis;

    catch
        error('One (or both) UV NO2 files missing')
    end

    gbs_no2uv=merge(ut,p);

    save([save_dir 'GBS_NO2_UV.mat'], 'gbs_no2uv')
end
end

function out = merge(t1_in,t2_in)

    % convert to arrays to allow manipulation
    t1=table2array(t1_in);
    t2=table2array(t2_in);

    % find matching twilights
    [~,ind1,ind2]=intersect(t1(:,1:3),t2(:,1:3),'rows');
    
    % average matching values (should have no NaNs in merged utgbs, pgbs files)
    out=(t1(ind1,:)+t2(ind2,:))/2;
    
    % replace error average with quadrature
    out(:,15)=sqrt( t1(ind1,15).^2 +t2(ind2,15).^2 )/2; % sigma_mean_vcd
    out(:,16)=sqrt( t1(ind1,16).^2 +t2(ind2,16).^2 )/2; % std_vcd
%     out(:,18)=sqrt( t1(ind1,18).^2 +t2(ind2,18).^2 )/2; % langley_vcd_err
    
    
    % add rest of data
    ind12 = setdiff(1:length(t1), ind1);
    ind22 = setdiff(1:length(t2), ind2);
    
    out=[out; t1(ind12,:); t2(ind22,:)];
    
    % sort by time
    out=sortrows(out,[1,2,3]);
    
    % convert back to table and add proper column names
    out=array2table(out);
    out.Properties.VariableNames=t1_in.Properties.VariableNames;
    
end



