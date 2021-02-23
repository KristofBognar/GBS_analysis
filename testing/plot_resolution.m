function plot_resolution()
% Plot resolutions from saved workspace file
% USAGE:
%   plot_resolution(inst)
% INPUT:
%   inst: 'P' for PEARL-GBS or 'U' for UT-GBS

% Testing folder needs to contain saved workspace (named resolution.mat)
% with variables:
%   {lamp}_{grating_center}_{var} where
% 
%       {lamp} = hg, ne or xe
%       {var}  = spec (contains loaded spectrum),
%                xpk  (contains x coord. of selected peaks), or
%                fwhm (contains FWHM of selected peaks)
%       {grating_center} = 1200_350, 600_450 or 300_450 for PEARL-GBS, or
%                          1800_350, 600_450 or 400_450 for UT-GBS
%      
% OUTPUT: Figures for resolution of all three gratings

load '/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/PEARL-GBSdata/PGBS_2020/testing_Eureka/res_tests/resolution.mat';


condition='P';

if condition=='P'
    
    %2019: exclude peaks not in last year's plots (neon peaks are problematic)
%     ne_600_450_xpk(1:2)=[];
%     ne_600_450_fwhm(1:2)=[];
%      ne_1200_350_xpk(4:5)=[];
%      ne_1200_350_fwhm(4:5)=[];
    
    % Plot 1200 grating
%     figure(1)
%     plot(hg_1200_350_xpk, hg_1200_350_fwhm.*0.043,'bx'),hold on
%     plot(ne_1200_350_xpk, ne_1200_350_fwhm.*0.043,'kx'),hold on
%     plot(xe_1200_350_xpk, xe_1200_350_fwhm.*0.043,'mx'),hold on
%     grid on
%     
%     legend('Hg', 'Ne', 'Xe')
%     xlim([0,2048])
%     ylim([0,1])
%     xlabel('CCD Pixels')
%     ylabel('FWHM (nm)')
% 
%     % Plot 600 grating
%     figure(2)
%     plot(hg_600_450_xpk, hg_600_450_fwhm.*0.107,'bx'),hold on
%     plot(ne_600_450_xpk, ne_600_450_fwhm.*0.107,'kx'),hold on
%     plot(xe_600_450_xpk, xe_600_450_fwhm.*0.107,'mx'),hold on
%     grid on
%     
%     legend('Hg', 'Ne', 'Xe')
%     xlim([0,2048])
%     ylim([0,3])
%     xlabel('CCD Pixels')
%     ylabel('FWHM (nm)')
% 
%     % Plot 300 grating
%     figure(3)
%     plot(hg_300_450_xpk, hg_300_450_fwhm.*0.208,'bx'),hold on
%     %plot(ne_300_450_xpk, ne_300_450_fwhm.*0.208,'kx'),hold on
%     plot(xe_300_450_xpk, xe_300_450_fwhm.*0.208,'mx'),hold on
%     grid on
% 
%     legend('Hg', 'Xe')
%     xlim([0,2048])
%     ylim([0,4])
%     xlabel('CCD Pixels')
%     ylabel('FWHM (nm)')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all on one figure, for masters report  
% 
    figure('Position', [100, 100, 680, 510])
    plot(hg_300_450_xpk, hg_300_450_fwhm.*0.208,'g+','linewidth',1.5,'DisplayName','300 gr/mm at 450 mm (2020)'),hold on
    plot(hg_600_450_xpk, hg_600_450_fwhm.*0.107,'ro','linewidth',1.5,'DisplayName','600 gr/mm at 450 nm (2020)'),hold on   
    plot(hg_1200_350_xpk, hg_1200_350_fwhm.*0.043,'bx','linewidth',1.5,'DisplayName','1200 gr/mm at 350 nm (2020)'),hold on
    
    plot(ne_1200_350_xpk, ne_1200_350_fwhm.*0.043,'bx','linewidth',1.5, 'HandleVisibility','off'),hold on
    plot(xe_1200_350_xpk, xe_1200_350_fwhm.*0.043,'bx','linewidth',1.5, 'HandleVisibility','off'),hold on
    
    %plot(ne_600_450_xpk(1), ne_600_450_fwhm(1)*0.107,'ro','linewidth',1.5),hold on
    %plot(ne_600_450_xpk(4), ne_600_450_fwhm(4)*0.107,'ro','linewidth',1.5),hold on
    plot(ne_600_450_xpk, ne_600_450_fwhm*0.107,'ro','linewidth',1.5, 'HandleVisibility','off'),hold on
    plot(xe_600_450_xpk, xe_600_450_fwhm.*0.107,'ro','linewidth',1.5, 'HandleVisibility','off'),hold on

    %plot(ne_300_450_xpk, ne_300_450_fwhm.*0.208,'g+','linewidth',1.5),hold on
    plot(xe_300_450_xpk, xe_300_450_fwhm.*0.208,'g+','linewidth',1.5, 'HandleVisibility','off'),hold on
    load '/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/GBS_testing_history/pgbs_2018/resolution.mat';
    plot(hg_1200_350_xpk, hg_1200_350_fwhm.*0.043,'kx','linewidth',1.5,'DisplayName','1200 gr/mm at 350 nm (2018)'),hold on
    plot(ne_1200_350_xpk, ne_1200_350_fwhm.*0.043,'kx','linewidth',1.5, 'HandleVisibility','off'),hold on
    plot(xe_1200_350_xpk, xe_1200_350_fwhm.*0.043,'kx','linewidth',1.5, 'HandleVisibility','off'),hold on

    load '/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/GBS_testing_history/pgbs_2016/resolution.mat';
    plot(hg_1200_350_xpk, hg_1200_350_fwhm.*0.043,'mx','linewidth',1.5,'DisplayName','1200 gr/mm at 350 nm (2016)'),hold on
    plot(ne_1200_350_xpk, ne_1200_350_fwhm.*0.043,'mx','linewidth',1.5, 'HandleVisibility','off'),hold on
    plot(xe_1200_350_xpk, xe_1200_350_fwhm.*0.043,'mx','linewidth',1.5, 'HandleVisibility','off'),hold on
    
    legend()


    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    

elseif condition=='U'
    
%     Plot 1800 grating
%     figure(1)
%     plot(hg_1800_350_xpk, hg_1800_350_fwhm.*0.03,'bx'),hold on
%     plot(ne_1800_350_xpk, ne_1800_350_fwhm.*0.03,'kx'),hold on
%     grid on
%     legend('Hg', 'Ne', 'Xe')
%     xlim([0,2048])
%     ylim([0,3])
%     xlabel('CCD Pixels')
%     ylabel('FWHM (nm)')
%     
% %     Plot 600 grating
%     figure(2)
%     plot(hg_600_450_xpk, hg_600_450_fwhm.*0.12,'bx'),hold on
%     plot(ne_600_450_xpk, ne_600_450_fwhm.*0.12,'kx'),hold on
%     plot(xe_600_450_xpk, xe_600_450_fwhm.*0.12,'mx'),hold on
%     grid on
%     legend('Hg', 'Ne', 'Xe')
%     xlim([0,2048])
%     ylim([0,3])
%     xlabel('CCD Pixels')
%     ylabel('FWHM (nm)')
%     
% %     Plot 400 grating
%     figure(3)
%     plot(hg_400_450_xpk, hg_400_450_fwhm.*0.16,'bx'),hold on
%     plot(ne_400_450_xpk, ne_400_450_fwhm.*0.16,'kx'),hold on
%     plot(xe_400_450_xpk, xe_400_450_fwhm.*0.16,'mx'),hold on
%     legend('Hg', 'Ne', 'Xe')
%     grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all on one figure, for masters report  
% 
%     figure('Position', [100, 100, 680, 510])
%     plot(hg_400_450_xpk, hg_400_450_fwhm.*0.16,'g+','linewidth',1.5,'DisplayName','400 gr/nm at 450 nm (2020)'),hold on
%     plot(hg_600_450_xpk, hg_600_450_fwhm.*0.12,'ro','linewidth',1.5,'DisplayName','600 gr/nm at 450 nm (2020)'),hold on
%     plot(hg_1800_350_xpk, hg_1800_350_fwhm.*0.03,'bx','linewidth',1.5,'DisplayName','1800 gr/nm at 350 nm (2020)'),hold on
%     
%     
%     plot(ne_1800_350_xpk, ne_1800_350_fwhm.*0.03,'bx','linewidth',1.5,'HandleVisibility','Off'),hold on
%    
%     plot(ne_600_450_xpk, ne_600_450_fwhm.*0.12,'ro','linewidth',1.5,'HandleVisibility','Off'),hold on
%     plot(xe_600_450_xpk, xe_600_450_fwhm.*0.12,'ro','linewidth',1.5,'HandleVisibility','Off'),hold on
% 
%     plot(ne_400_450_xpk, ne_400_450_fwhm.*0.16,'g+','linewidth',1.5,'HandleVisibility','Off'),hold on
%     plot(xe_400_450_xpk, xe_400_450_fwhm.*0.16,'g+','linewidth',1.5,'HandleVisibility','Off'),hold on
%     load '/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/GBS_testing_history/utgbs_2018/resolution.mat';
%     plot(hg_600_450_xpk, hg_600_450_fwhm.*0.12,'ko','linewidth',1.5,'DisplayName','600 gr/nm at 450 nm (2018)'),hold on
%     plot(ne_600_450_xpk, ne_600_450_fwhm.*0.12,'ko','linewidth',1.5,'HandleVisibility','Off'),hold on
%     plot(xe_600_450_xpk, xe_600_450_fwhm.*0.12,'ko','linewidth',1.5,'HandleVisibility','Off'),hold on
%     load '/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/GBS_testing_history/utgbs_2016/resolution.mat';
%     plot(hg_600_450_xpk, hg_600_450_fwhm.*0.12,'mo','linewidth',1.5,'DisplayName','600 gr/nm at 450 nm (2016)'),hold on
%     plot(ne_600_450_xpk, ne_600_450_fwhm.*0.12,'mo','linewidth',1.5,'HandleVisibility','Off'),hold on
%     plot(xe_600_450_xpk, xe_600_450_fwhm.*0.12,'mo','linewidth',1.5,'HandleVisibility','Off'),hold on
%     
%     legend()
%       
    

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    xlim([0,2048])
    ylim([0,3])
    xlabel('CCD Pixels')
    ylabel('FWHM (nm)')
    
else
    disp('---')
    disp('Error: Invalid instrument (see documentation below)')
    disp('---')
    help plot_resolution

end


    
    

