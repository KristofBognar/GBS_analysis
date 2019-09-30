function plot_oclo_dscds(year)
% Plot and save OClO dscds for given year
%
% Requires:
%    QDOAS output in <year>/output/oclo_dscd.ASC 

%%

if year<=2011, 
    source='C'; 
else
    source='asd';
end

% change to working directory
cd /home/kristof/work/GBS/PEARL-GBS;

load 2016/avg_twilight_for_oclo/vcd_input.mat

file_nm=[num2str(year),'/output/oclo_dscd.ASC'];
if source=='C', file_nm=['pre-2011_oclo/oclo_cristen_',num2str(year),'.ASC']; end

[dscd_S,qdoas_filt,qdoas_raw]=read_QDOAS_oclo(file_nm,col_no2_p0,filt_good,1,source);

% filter out unvanted SZA
ind_good=find(dscd_S.sza>=86 & dscd_S.sza<=91);

% filter out cold weekend for 2016
if year==2016
    ind_good=find((dscd_S.sza>=86 & dscd_S.sza<=91) & (dscd_S.fd<72 | dscd_S.fd>78.74));
end

qdoas_filt2=qdoas_filt(ind_good,:);


fd=[];
oclo_dscd=[];
oclo_dscd_std=[]; % random error
oclo_err=[]; % systematic error
oclo_rms=[]; % systematic error

% loop ower days
for i=floor(dscd_S.fd(1)):ceil(dscd_S.fd(end))
    % loop over both twilights
    for j=1:2

        % noon at .74 fd
        % midnight at .24
        if j==1
            lower=i+0.24;
            upper=i+0.74;

            % find spectra for twilight
            ind_am=find(qdoas_filt2(:,col_no2_p0.sza)>=89 & qdoas_filt2(:,col_no2_p0.sza)<=91 & ...
                        qdoas_filt2(:,col_no2_p0.fd)>=lower & qdoas_filt2(:,col_no2_p0.fd)<upper);
            
        elseif j==2
            lower=i+0.74;
            upper=i+1.24;

            % find spectra for twilight
            ind_pm=find(qdoas_filt2(:,col_no2_p0.sza)>=89 & qdoas_filt2(:,col_no2_p0.sza)<=91 & ...
                        qdoas_filt2(:,col_no2_p0.fd)>=lower & qdoas_filt2(:,col_no2_p0.fd)<upper);
            
        end
    end
             
    % assign results to arrays (1st row am, 2nd row pm)    
    if ~isempty(ind_am) || ~isempty(ind_pm)
       fd=[fd, [nanmean(qdoas_filt2(ind_am,col_no2_p0.fd));...
                nanmean(qdoas_filt2(ind_pm,col_no2_p0.fd))]]; 
       oclo_dscd=[oclo_dscd, [nanmean(qdoas_filt2(ind_am,col_no2_p0.oclo_dscd)); ...
                  nanmean(qdoas_filt2(ind_pm,col_no2_p0.oclo_dscd))]]; 
       oclo_err=[oclo_err, [nanmean(qdoas_filt2(ind_am,col_no2_p0.oclo_err)); ...
                  nanmean(qdoas_filt2(ind_pm,col_no2_p0.oclo_err))]]; 
       oclo_rms=[oclo_rms, [nanmean(qdoas_filt2(ind_am,col_no2_p0.rms)); ...
                  nanmean(qdoas_filt2(ind_pm,col_no2_p0.rms))]]; 
       oclo_dscd_std=[oclo_dscd_std, [nanstd(qdoas_filt2(ind_am,col_no2_p0.oclo_dscd)); ...
                  nanstd(qdoas_filt2(ind_pm,col_no2_p0.oclo_dscd))]]; 
    end
end

% calculate total error (quadrature of syst. and random errors)
error=(oclo_dscd_std.^2 + oclo_err.^2 + oclo_rms.^2).^0.5;
[fd_tot,ind]=sort([fd(1,:),fd(2,:)]);
error=[error(1,:),error(2,:)];
error=error(ind);

figure(10)
hold on
grid on

plot(fd(1,:),oclo_dscd(1,:),'b.', 'markersize',16)
plot(fd(2,:),oclo_dscd(2,:),'r.', 'markersize',16)

% plot detection limit (3*sigma)
% plot(fd_tot, error.*3, 'k--')
det_lim=nanmean(error)*3;
plot([40,125],[det_lim,det_lim], 'k-')

% xlim([65,110])

legend('a.m.', 'p.m.', 'Approx. detection limit')

xlabelstr=['Day of the Year, ', num2str(year),' (UTC)'];
xlabel(xlabelstr)
ylabel('OClO DSCD (molec/cm^2)')
% title('OClO DSCDs (89<=SZA<91) from Eureka, March 5 - April 20')

% figure(2)
% hold on
% 
% plot(fd(1,:),oclo_dscd_std(1,:),'b.-')
% plot(fd(2,:),oclo_dscd_std(2,:),'r.-')
% 
% plot(fd(1,:),oclo_err(1,:),'b.--')
% plot(fd(2,:),oclo_err(2,:),'r.--')
% 
% plot(fd_tot,error,'k--')

%% write results to file

% get total integration time (in seconds)
tint=qdoas_filt2(:,5).*qdoas_filt2(:,6);
%create array for viewing elevation
elev=ones(size(dscd_S.fd(ind_good))).*90.0;
% convert SAA to north=0, east=90
dscd_S.saa=dscd_S.saa+180;

% create output array
temp=double(struct2dataset(dscd_S));
% our fd starts with 1 on jan. 1, 00:00
out_array=[dscd_S.year(ind_good),dscd_S.fd(ind_good)-1,tint,elev,temp(ind_good,4:end-1)];

% write header in the file
f_out=[num2str(year),'/avg_twilight_for_oclo/OClO_Eureka_',num2str(year),'_UToronto.asc'];
if source=='C', f_out=['pre-2011_oclo/OClO_Eureka_',num2str(year),'_UToronto.asc']; end
fid = fopen(f_out, 'w');

fprintf(fid(1), '%s\n', 'PEARL Ground Based Spectrometer (PEARL-GBS) zenith-sky data');
fprintf(fid(1), '%s\n', 'Location: Eureka, Nunavut, Canada (80.053N, 86.416W)');
fprintf(fid(1), '%s\n', 'PI: Kimberly Strong (email: strong@atmosp.physics.utoronto.ca)');
fprintf(fid(1), '%s\n', 'University of Toronto, 60 St George St., Toronto, Ontario, M5S 1A7, Canada');
fprintf(fid(1), '%s\n', 'File created by: Kristof Bognar (email: kbognar@physics.utoronto.ca)');
fprintf(fid(1), '%s\n', '');

if source=='C'
    fprintf(fid(1), '%s\n', 'Retrieval done by a former graduate student (see Adams et al., 2012)');
    fprintf(fid(1), '%s\n', 'Fitting Window: 350-380 nm');
    fprintf(fid(1), '%s\n', 'Polynomial: 3rd order');
    fprintf(fid(1), '%s\n', 'Offset: 1st order');
    fprintf(fid(1), '%s\n', 'Cross-sections: OClO (204K, Wahner et al., 1987), NO2 (220 K, Vandaele et al., 1997), O3 (223 K, Bugomil et al., 2003), O4 (Greenblatt et al., 1990), BrO (223K, Fleischmann et al., 2004), Ring (Chance and Spurr, 1997)');
    fprintf(fid(1), '%s\n', '*');
    fprintf(fid(1), '%s\n', 'Twilight spectra were averaged to a time-resolution of 20-30 minutes  prior to QDOAS analysis');
    fprintf(fid(1), '%s\n', '');
else
    fprintf(fid(1), '%s\n', 'Retrieval code: QDOAS V2.109.4');
    fprintf(fid(1), '%s\n', 'Fitting Window: 350-380 nm');
    fprintf(fid(1), '%s\n', 'Polynomial: 3rd order');
    fprintf(fid(1), '%s\n', 'Offset: 1st order');
    fprintf(fid(1), '%s\n', 'Cross-sections: OClO (204K, Wahner et al., 1987), NO2 (220 K, Vandaele et al., 1997), O3 (223 K, Bugomil et al., 2003), O4 (Hermans et al., 2003), BrO (223K, Fleischmann et al., 2004), Ring (Chance and Spurr, 1997)');
    fprintf(fid(1), '%s\n', 'Reference: noon zenith spectrum averaged in a 0.03 degree SZA window');
    fprintf(fid(1), '%s\n', 'Twilight spectra were averaged in 0.5 degree SZA bins prior to QDOAS analysis');
    fprintf(fid(1), '%s\n', '');
end

head=sprintf('Year\tDOY\tTint(s)\tViewing_Elev\tSZA\tSAA\tRMS\tOClO_DSCD\tOClO_DSCD_error\tSpectrum_shift\tSpectrum_stretch\tReference_SZA');
fprintf(fid(1), '%s\n', head);

fclose(fid);

% Writre data in the file
dlmwrite(f_out,out_array,'delimiter','\t','precision',8,'-append')




