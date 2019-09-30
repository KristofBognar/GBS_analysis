% function plot_QDOAS_rms
% To plot QDOAS rms errors. Has to run from VCD directory, must load QDOAS
% files manually

% check if qdoas results are loaded
if ~exist('qdoas_raw'), error('Must load QDOAS file to be plotted'), return, end

% load filtering and column info
load('../vcd_input.mat')


figure(1)

if ~isempty(strfind(pwd,'UT-GBS'))
    
    subplot(211)

%     ind=find(qdoas_raw(:,col_o3_3.sza) <=90 & qdoas_raw(:,col_o3_3.sza) >= 86);
    ind=find(qdoas_raw(:,col_o3_3.sza) <=90);
    
    plot(qdoas_raw(ind,col_o3_3.fd), qdoas_raw(ind,col_o3_3.rms), 'k.'), hold on
%     xlim([54,last_day])
    ylim([0,0.005])
    ylabel('Ozone RMS')
    title('UT-GBS, SZA < 90')

    subplot(212)

%     ind=find(qdoas_raw(:,col_no2_3.sza) <=91 & qdoas_raw(:,col_no2_3.sza) >= 86);
    ind=find(qdoas_raw(:,col_no2_3.sza) <=91);

    plot(qdoas_raw(ind,col_no2_3.fd), qdoas_raw(ind,col_no2_3.rms), 'k.'), hold on
%     xlim([54,last_day])
    ylim([0,0.005])
    xlabel('Freactional day (UTC)')
    ylabel('NO2 RMS')
    title('UT-GBS, SZA < 91')
    
elseif ~isempty(strfind(pwd,'PEARL-GBS'))
%     subplot(211)
%     ind=find(qdoas_raw(:,col_no2_p0.sza) <=91 & qdoas_raw(:,col_no2_p0.sza) >= 86);
    ind=find(qdoas_raw(:,col_no2_p0.sza) <=91);
%     [maxima,maxind] = findpeaks(qdoas_raw(:,col_no2_p0.sza));

    plot(qdoas_raw(ind,col_no2_p0.fd), qdoas_raw(ind,col_no2_p0.rms), 'k.'), hold on
%     xlim([54,last_day])
    ylim([0,0.005])
    xlabel('Freactional day (UTC)')
    ylabel('NO2 RMS')
    title('PEARL-GBS, SZA < 91')
    
%     subplot(212)
%     plot(qdoas_raw(:,col_no2_p0.fd), qdoas_raw(:,col_no2_p0.sza), 'k-'), hold on


else
    disp('Must specify instrument: P for PEARL-GBS, U for UT-GBS');
end