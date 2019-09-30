

instr='UT-GBS';
% instr='PEARL-GBS';

% compare='dscd';
compare='vcd';


% for year=2000:2000
year=2003;
% if year==2012, continue, end

if strcmp(compare,'dscd')

    
    
    if strcmp(instr,'UT-GBS') && year<2017
        merged='_merged';
    else
        merged='';
    end
     
    % load QDOAS table
    tablename=['/home/kristof/work/GBS/QDOAS_results/yearly_tables/' ...
               instr '_' num2str(year) merged '.mat'];
    try
        load(tablename);
    catch
%         continue
    end

    ind=find(data.NO2SlColno2<1e19 & data.NO2_425450SlColno2<1e19);
%     ind=find(data.NO2SlColno2<1e19 & data.NO2_425450SlColno2<1e19 & data.SZA<91);
    
% %     figure(1)
% %     ax1=subplot(311);
% %     plot(data.Fractionalday(ind),data.NO2_425450SlColno2(ind),'ro'), hold on
% %     plot(data.Fractionalday(ind),data.NO2SlColno2(ind),'kx'), hold on
% %     ylabel('NO_2 dSCD')
% %     
% %     ax2=subplot(312);
% %     plot(data.Fractionalday(ind),data.NO2_425450RMS(ind),'ro'), hold on
% %     plot(data.Fractionalday(ind),data.NO2RMS(ind),'kx'), hold on
% %     ylabel('NO_2 RMS')
% %     
% %     ax3=subplot(313);
% %     plot(data.Fractionalday(ind),data.NO2_425450SlErrno2(ind),'ro'), hold on
% %     plot(data.Fractionalday(ind),data.NO2SlErrno2(ind),'kx'), hold on
% %     ylabel('NO_2 fit error')
% %     xlabel(['Day of the year, ' num2str(year)])
% %     
% %     hlink = linkprop([ax1,ax2,ax3],{'xlim'});     
    
    dscd_diff=mean_diff('rel',data.NO2SlColno2(ind),data.NO2_425450SlColno2(ind));
    rms_diff=mean_diff('rel',data.NO2RMS(ind),data.NO2_425450RMS(ind));
    err_diff=mean_diff('rel',data.NO2SlErrno2(ind),data.NO2_425450SlErrno2(ind));
    
% %     disp(num2str(year))
% %     disp(['dscd: ' num2str(dscd_diff) '    rms: ' num2str(rms_diff) ...
% %           '    err: ' num2str(err_diff)])
% %     
% %     pause
% %     
% %     clf
    
    figure(2)
    
    subplot(311)
    plot(data.Year(1),dscd_diff,'kx'), hold on
    ylabel('\Delta dSCD')
    subplot(312)
    plot(data.Year(1),rms_diff,'kx'), hold on
    ylabel('\Delta RMS')
    subplot(313)
    plot(data.Year(1),err_diff,'kx'), hold on
    ylabel('\Delta fit error')
    xlabel('Year')

    clearvars data ind
    

elseif strcmp(compare,'vcd')
    
    % load small window VCDs
    load('/home/kristof/work/GBS/VCD_results/UT-GBS_reanalysis_old_err_budget/UT-GBS_NO2_VCD_all.mat');
    
% %     % load new VCDs table
% %     filename=['/home/kristof/work/GBS/VCD_results/' instr '_NO2_VCD_' num2str(year) '.mat'];
% %     try
% %         load(filename);
% %     catch
% %         continue
% %     end
% %     
% %     [~,ind_small,ind_large] = intersect(reanalysis(:,1:3),VCD_table(:,1:3));
% %     
% %     vcd_diff=mean_diff('rel',VCD_table.mean_vcd(ind_large),reanalysis.mean_vcd(ind_small));
% % 
% % % %     figure(1)
% % % %     plot(VCD_table.year(1),vcd_diff,'kx'), hold on
% % % %     ylabel('\Delta VCD')
% %     
% %     figure(2)
% %     gscatter(reanalysis.fd(ind_small), reanalysis.mean_vcd(ind_small), ...
% %              reanalysis.ampm(ind_small), 'rb', '..'), hold on
% %     gscatter(VCD_table.fd, VCD_table.mean_vcd, ...
% %              VCD_table.ampm, 'rb', 'xx'), hold off

    figure(2)
    ind_small=find(reanalysis.year==year);
    gscatter(reanalysis.fd(ind_small), reanalysis.mean_vcd(ind_small), ...
             reanalysis.ampm(ind_small), 'rb', '..'), hold on
         
    % load large window VCDs
    load('/home/kristof/work/GBS/VCD_results/UT-GBS_NO2_VCD_all.mat');
    ind_large=find(reanalysis.year==year);
    gscatter(reanalysis.fd(ind_large), reanalysis.mean_vcd(ind_large), ...
             reanalysis.ampm(ind_large), 'rb', 'xx'), hold off

    
    
end
% end
