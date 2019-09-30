% plot RMS and DSCDs for fitting windows with different settings

%% input
file_nm='/home/kristof/work/GBS/UT-GBS/2011/output/UT-GBS_2011_reanalysis.csv';
col=create_colnum_struct('gbs_redo');


%% Read in the QDOAS data
qdoas_raw = csvread(file_nm);

% sort by fractional day and ID whether there are doubles of some values
qdoas_raw = sortrows(qdoas_raw, col.fd);
all_ind = 1:length(qdoas_raw(:,1));
[a, unique_ind,b] = unique(qdoas_raw(:, col.fd));
diff_ind = setdiff(all_ind, unique_ind);
if ~isempty(diff_ind)
    disp('[WARNING]: File contains multiple entries taken at the same time.')
    disp('frac_day')
    disp('---------')
    disp(qdoas_raw(diff_ind, col.fd))
end

% % filter out rms = 0, and dscd = 9999
% ind = find(qdoas_raw(:, col.rms) ~= 0 & ...
%     qdoas_raw(:, col.dscd) ~= 9999 & ...
%     qdoas_raw(:, col.err) ~= 9999 & ...
%     qdoas_raw(:, col.sza) < filt.sza_max);
% qdoas = qdoas_raw(ind,:);


%% plot stuff

fd=qdoas_raw(:,col.fd);

% compare stuff
% ind=find(qdoas_raw(:,col.o3_rms) <9999 & qdoas_raw(:,col.o3_offline_rms) <9999 &...
%          qdoas_raw(:,col.o3_rms) >0 & qdoas_raw(:,col.o3_offline_rms) >0 &...
%          qdoas_raw(:,col.o3_293_rms) >0 & qdoas_raw(:,col.o3_293_rms) <9999 &...
%          qdoas_raw(:,col.sza) < 91);
     
ind=find(qdoas_raw(:,col.no2_rms) <9999 & qdoas_raw(:,col.no2_298_rms) <9999 &...
         qdoas_raw(:,col.no2_rms) >0 & qdoas_raw(:,col.no2_298_rms) >0 &...
         qdoas_raw(:,col.sza) < 91);

w_c_o3=(qdoas_raw(ind,col.o3_293_rms)./qdoas_raw(ind,col.o3_rms)-1)*100;
w_c_no2=(qdoas_raw(ind,col.no2_298_dscd)./qdoas_raw(ind,col.no2_dscd)-1)*100;
on_off=(qdoas_raw(ind,col.o3_offline_rms)./qdoas_raw(ind,col.o3_rms)-1)*100;
     
plot(fd(ind),w_c_no2,'k.')



