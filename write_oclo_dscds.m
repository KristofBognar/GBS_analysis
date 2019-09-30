

% code to format 2017 and later OClO dSCDs from QDOAS
% file format based on what I sent to Gaia in 2017 (files were probably
% straight from QDOAS)

% save dolder
savedir='/home/kristof/work/documents/GBS_data_for_people/Gaia_OClO/';


% load data (saved manually)
load('/home/kristof/work/GBS/QDOAS_results/oclo_dscds_2017-2019.mat');


for yr=2017:2019
    
    %% data processing
    % get correct variable
    data=eval(['oclo_' num2str(yr)]);
    
    % filter data (use NO2-UV, RMS thresholds of 0.002 for SZA<87 and 0.003 for SZA>=87)
    ind_filt=find((data.NO2RMS>=0.002 & data.SZA<87) |...
                  (data.NO2RMS>=0.003 & data.SZA>=87));
              
    data(ind_filt,:)=[];
    data(data.SZA>92 | data.SZA<86,:)=[];
    
    
    % construct output array
    out_array=NaN(size(data,1),12);
    
    out_array(:,1)=data.Year;
    out_array(:,2)=data.Fractionalday-1; % qdoas output is 1 on jan 1, 00:00
    out_array(:,3)=data.Tint;
    out_array(:,4)=data.Elevviewingangle;
    out_array(:,5)=data.SZA;
    out_array(:,6)=data.SolarAzimuthAngle;
    out_array(:,7)=data.NO2RMS;
    out_array(:,8)=data.NO2SlColoclo;
    out_array(:,9)=data.NO2SlErroclo;
    out_array(:,10)=data.NO2ShiftSpectrum;
    out_array(:,11)=data.NO2StretchSpectrum1;
    out_array(:,12)=data.NO2RefZm;
    
    %% write output
    % output file
    f_out=[savedir 'OCLO_Eureka_' num2str(yr) '_UToronto.asc'];
    
    % open file
    fid = fopen(f_out, 'w');

    % write header
    fprintf(fid(1), '%s\n', '*PEARL Ground Based Spectrometer (PEARL-GBS) zenith-sky data');
    fprintf(fid(1), '%s\n', '*Location: Eureka, Nunavut, Canada (80.053N, 86.416W)');
    fprintf(fid(1), '%s\n', '*PI: Kimberly Strong (email: strong@atmosp.physics.utoronto.ca)');
    fprintf(fid(1), '%s\n', '*University of Toronto, 60 St George St., Toronto, Ontario, M5S 1A7, Canada');
    fprintf(fid(1), '%s\n', '*File created by: Kristof Bognar (email: kbognar@physics.utoronto.ca)');
    fprintf(fid(1), '%s\n', '*');
    fprintf(fid(1), '%s\n', '*Retrieval code: QDOAS V3.1');
    fprintf(fid(1), '%s\n', '*Fitting Window: 350-380 nm');
    fprintf(fid(1), '%s\n', '*Polynomial: 3rd order');
    fprintf(fid(1), '%s\n', '*Offset: 1st order');
    fprintf(fid(1), '%s\n', '*Cross-sections: OClO (204K, Wahner et al., 1987), NO2 (220 K, Vandaele et al., 1997), O3 (223 K, Bugomil et al., 2003), O4 (Greenblatt et al., 1990), BrO (223K, Fleischmann et al., 2004), Ring (Chance and Spurr, 1997)');
    fprintf(fid(1), '%s\n', '*Reference: daily noon zenith spectrum averaged in a 0.03 degree SZA window');
    fprintf(fid(1), '%s\n', '*Twilight spectra were averaged in 0.5 degree SZA bins prior to QDOAS analysis');
    fprintf(fid(1), '%s\n', '*');
    fprintf(fid(1), '%s\n', '*DOY (UTC): January 1st, 00:00 = 0');
    fprintf(fid(1), '%s\n', '*');
    
    % header
    head=sprintf('*Year\tDOY\tTint(s)\tViewing_Elev\tSZA\tSAA\tRMS\tOClO_DSCD\tOClO_DSCD_error\tSpectrum_shift\tSpectrum_stretch\tReference_SZA');
    fprintf(fid(1), '%s\n', head);
    
    fclose(fid);
    
    % data
    dlmwrite(f_out,out_array,'delimiter','\t','precision',8,'-append')

    % plot results
    errorbar(data.Fractionalday,data.NO2SlColoclo,data.NO2SlErroclo,...
             'linestyle','none'), hold on
         
end







% NO 2 -vis, RMS thresholds of 0.002 for SZA < 90° and 0.003 for SZA > 90° were
% selected to remove most extreme outlying DSCDs and DSCDs retrieved from spectra that were
% visually identified as saturated. For SZA < 94°, more than 95% of the data passed this RMS
% threshold, as shown in Table 4.2. For NO 2 -UV, weaker RMS thresholds of 0.002 for SZA < 87°
% and 0.003 for SZA < 91° were used