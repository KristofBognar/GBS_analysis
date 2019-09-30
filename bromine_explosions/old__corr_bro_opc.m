% correlation plots for BrO VCDs / AOD and OPC supermicron particle count
% old code, newer files exist that do similar things

%% final variables

vcd_aod_all=[];
vcd_top_all=[];
ft_prof_all=[];
ft_opc_all=[];
sm_all=[];
subm_all=[];
wspd_all=[];
T_inv=[];


%% loop over years
for year=2016:2017;
    %% load BrO/AOD data
    
    % select aerosol or tracegas
    option='tg';
    
    % select layer to stat partial VCD/AOD from
    top=5;
    
    if year==2016
        %59 for leap years (2016), since data is already fractional day
        subtract=59; 
    else
        subtract=58;
    end

    if option=='a'
        filedir=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
                 num2str(year) '/aerosol/'];
    elseif option=='tg'
        filedir=['/home/kristof/work/profile_retrievals/profile_results/eureka_'...
                 num2str(year) '/tracegas/'];
    end

    % list of matlab files with daily profile data
    if year==2017
        files={'20170307','20170308','20170309','20170310','20170311','20170312',...
               '20170313','20170314','20170315','20170316','20170317','20170318',...
               '20170319','20170320'};
    elseif year==2016 && strcmp(option,'tg')
        files={'20160319_noadaptive_1iter_skip0.mat',...
               '20160320_noadaptive_1iter_skip0.mat',...
               '20160321_noadaptive_1iter_skip0.mat',...
               '20160322_noadaptive_1iter_skip0.mat'};
    elseif year==2016 && strcmp(option,'a')
        files={'20160319_noadaptive_skip0.mat',...
               '20160320_noadaptive_skip0.mat',...
               '20160321_noadaptive_skip0.mat',...
               '20160322_noadaptive_skip0.mat'};
    end
    % files={'20160319_noadaptive-1.mat',...
    %        '20160320_noadaptive-1.mat',...
    %        '20160321_noadaptive-1.mat',...
    %        '20160322_noadaptive-1.mat'};

    % files={'20160319_noadaptive_no0.mat',...
    %        '20160320_noadaptive_no0.mat',...
    %        '20160321_noadaptive_no0.mat',...
    %        '20160322_noadaptive_no0.mat'};

    % files={'20160319_noadaptive_no0_1iter.mat',...
    %        '20160320_noadaptive_no0_1iter.mat',...
    %        '20160321_noadaptive_no0_1iter.mat',...
    %        '20160322_noadaptive_no0_1iter.mat'};

    % variable columns (files are per day):
    %   info: DoF, VCD/extinction, error, fractional time
    %   prof: rows are altitude, columns are time, units are ppm
    %   prof_nd: rows are altitude, columns are time, units are molec/cm^3
    %   prof_err: same as profiles
    %   dscd: SZA,  elev,  rel_azim,  wl,  BrOmeas,  err_BrOmeas,  BrOretr (or O4 for aerosol)

    % limit=0.25; % for aerosol extinction

    % variables for VCD/AOD and time
    vcd_aod=[];
    vcd_top=[];
    ft_prof=[];

    for i=1:size(files,2)

        % load file 
        load([filedir,files{i}])

        if option=='a'

            % replace high values
            if exist('limit','var')
                prof(prof>limit)=limit;
            end
            % only plot good profiles, filter by DoF
    %         ind=find(info(:,1)>1.5);
            ind=find(info(:,1)>0);
            
            part_prof=prof*0.2; % layers are 0.2 km thick

        elseif option=='tg'

            % convert to ppt (results are in ppm)
            prof=prof*1e6;
            % only plot good profiles, filter by DoF
    %         ind=find(info(:,1)>0.8);
            ind=find(info(:,1)>0);

            part_prof=prof_nd(:,ind)*20000; % layers are 200m thick
    %         part_prof_err=prof_nd_err(:,ind)*20000;
        end
        
        % calculate partial profile aloft
        for jj=1:size(prof,2)
            vcd_top=[vcd_top; sum(part_prof(top:end,jj))];
        end

        % asign variables
        vcd_aod=[vcd_aod; info(:,2)];
        ft_prof=[ft_prof; info(:,4)];

    end

    %% load OPC data

    if year==2016
        load(['/home/kristof/work/SMPS/opc_size_dist_0316.mat']);
    elseif year==2017
        load(['/home/kristof/work/SMPS/opc_size_dist_0317.mat']);
    end
    
    % convert times to fractional time
    [ft_opc]=fracdate(time,'dd-eee-yyyy hh:mm:ss');

    % bins to include
    b1=3;
    b2=5;

    % convert dN/dlog10(Dp) to N
    Dp=[Dp,20000];
    logDp=log10(Dp(2:end)./Dp(1:end-1));
    for i=1:length(logDp)
        Dp_data(:,i)=Dp_data(:,i)*logDp(i);
    end
    
    % get supermicron particle info
    supermicron_full=sum(Dp_data(:,b1:b2),2);
    supermicron=[];

    % filter out bad data
    ind=find(supermicron_full==0 | supermicron_full>400);
    supermicron_full(ind)=[];
    ft_opc(ind)=[];

    %% load windspeed data from PWS
    
    load(['/home/kristof/work/weather_stations/ridge_lab/PWS_03_' num2str(year) '.mat']);
    wspd_tmp=[];
    
    %% load and integrate smps data
    if year==2016
        load(['/home/kristof/work/SMPS/size_dist.mat']);
    elseif year==2017
        load(['/home/kristof/work/SMPS/size_dist_0317.mat']);
    end
    
    % convert times to fractional time
    [ft_smps]=fracdate(time,'dd-eee-yyyy hh:mm:ss');
    
    % convert dN/dlog10(Dp) to N
    Dp_edge=(Dp(2:end)+Dp(1:end-1))/2;
    Dp_edge=[10,Dp_edge,500];
    logDp=log10(Dp_edge(2:end)./Dp_edge(1:end-1));
    for i=1:length(logDp)
        Dp_data(:,i)=Dp_data(:,i)*logDp(i);
    end
    
    
    % integrate Dp data
    % 50-500 nm: 23-end
    % 100-500 nm: 33-end
    submicron_full=sum(Dp_data(:,33:end),2);
%     submicron_full=tot_data;
    submicron=[];
    
    
    %% average stuff +- 7 min around profile time
    % (about the length of each scan)
    % convert min to fractional time
    bin=7/(60*24);

    
    for i=1:length(ft_prof)
        
        % average particle concentration from OPC
        ind=find(ft_opc>ft_prof(i)-bin & ft_opc<ft_prof(i)+bin);
        supermicron=[supermicron; nanmean(supermicron_full(ind))];

        % average particle concentration from SMPS
        ind=find(ft_smps>ft_prof(i)-bin & ft_smps<ft_prof(i)+bin);
        submicron=[submicron; nanmean(submicron_full(ind))];

        % average windspeed from PWS
        ind=find(ft_wnd>ft_prof(i)-bin & ft_wnd<ft_prof(i)+bin);
        wspd_tmp=[wspd_tmp; nanmean(wspd(ind))];
        
        % temperature inversion data (visual inspection)
        if year==2016 && any([21,22]==floor(ft_prof(i)-59))
%             T_inv=[T_inv; 1]; % uncomment if you want t inversion info on plot
            T_inv=[T_inv; 0];
        elseif year==2016 && any([19,20]==floor(ft_prof(i)-59))
            T_inv=[T_inv; 0];
        elseif year==2017 && any([8,9,10,11,12,13,19,20]==floor(ft_prof(i)-58))
%             T_inv=[T_inv; 1];
            T_inv=[T_inv; 0];
        elseif year==2017 && any([7,14,15,16,17,18]==floor(ft_prof(i)-58))
            T_inv=[T_inv; 0];
        end
    end

    % do something about weird kink in OPC data for 2017
    if year==2017
        supermicron(195:197)=interp1([ft_prof(194),ft_prof(198)],...
                                     [supermicron(194),supermicron(198)],...
                                     ft_prof(195:197)).*[0.9;0.95;1.03];
    end
    
    % remove profiles with no OPC data
%     vcd_aod(ft_prof>ft_opc(end))=[];
%     vcd_top(ft_prof>ft_opc(end))=[];
%     ft_prof(ft_prof>ft_opc(end))=[];
    
    %% assign final variables
    vcd_aod_all=[vcd_aod_all; vcd_aod];
    vcd_top_all=[vcd_top_all; vcd_top];
    ft_prof_all=[ft_prof_all; ft_prof];
    ft_opc_all=[ft_opc_all; ft_opc];
    sm_all=[sm_all;supermicron];
    subm_all=[subm_all;submicron];
    wspd_all=[wspd_all;wspd_tmp];
    
end


% figure(1)                         
% subplot(211)
% plot(ft_prof, vcd_aod, 'k.'), hold on
% subplot(212)
% plot(ft_prof-59, supermicron, 'ko'), hold on
% plot(ft_opc-59, supermicron_full, 'b-'), hold on
% % xlim([65,80])
% xlim([19.5,23])
% ylim([0,500])

% figure(2)
% plot(sm_all, vcd_aod_all, 'k.')
% 

ind=~isnan(wspd_all);
T_inv=logical(T_inv);

[~,tmp]=min(wspd_all);
wspd_all(tmp)=0.01;

wspd_all(wspd_all>14)=14;


% Create figure
figure1 = figure;
colormap(jet(14));

%% plot upper half of VCD vs supermicrom particle conc.
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.6 0.33 0.3]);
hold(axes1,'on');

% scatter(sm_all(~ind),vcd_top_all(~ind),25,'k','filled');
% scatter(sm_all(ind),vcd_top_all(ind),25,wspd_all(ind),'filled');

% with temperature inversion below PEARL
scatter(sm_all(~ind & T_inv),vcd_top_all(~ind & T_inv),50,'k','filled','marker','pentagram');
scatter(sm_all(ind & T_inv),vcd_top_all(ind & T_inv),50,wspd_all(ind & T_inv),'filled','marker','pentagram');

% no temperature inversion below PEARL
scatter(sm_all(~ind & ~T_inv),vcd_top_all(~ind & ~T_inv),25,'k','filled');
scatter(sm_all(ind & ~T_inv),vcd_top_all(ind & ~T_inv),25,wspd_all(ind & ~T_inv),'filled');

ylim([0,4.5e13])

box on
% Create ylabel
ylabel('BrO VCD_{0.7-4km} (mol/cm^2)');

%% plot upper half of VCD vs submicron particle conc.
axes2 = axes('Parent',figure1,...
    'Position',[0.55 0.6 0.33 0.3]);
hold(axes2,'on');

scatter(subm_all(~ind & T_inv),vcd_top_all(~ind & T_inv),50,'k','filled','marker','pentagram');
scatter(subm_all(ind & T_inv),vcd_top_all(ind & T_inv),50,wspd_all(ind & T_inv),'filled','marker','pentagram');

scatter(subm_all(~ind & ~T_inv),vcd_top_all(~ind & ~T_inv),25,'k','filled');
scatter(subm_all(ind & ~T_inv),vcd_top_all(ind & ~T_inv),25,wspd_all(ind & ~T_inv),'filled');

ylim([0,4.5e13])
% xlim([2,12]*1000)

box on

%% plot full VCD vs supermicrom particle conc.
axes3 = axes('Parent',figure1,...
    'Position',[0.13 0.22 0.33 0.3]);
hold(axes3,'on');

scatter(sm_all(~ind & T_inv),vcd_aod_all(~ind & T_inv),50,'k','filled','marker','pentagram');
scatter(sm_all(ind & T_inv),vcd_aod_all(ind & T_inv),50,wspd_all(ind & T_inv),'filled','marker','pentagram');

scatter(sm_all(~ind & ~T_inv),vcd_aod_all(~ind & ~T_inv),25,'k','filled');
scatter(sm_all(ind & ~T_inv),vcd_aod_all(ind & ~T_inv),25,wspd_all(ind & ~T_inv),'filled');

ylim([0,12e13])

box on
% Create xlabel
xlabel('D_p > 1\mum (cm^-^3)');

% Create ylabel
ylabel('BrO VCD_{0-4km} (mol/cm^2)');

%% plot full VCD submicron particle conc.
axes4 = axes('Parent',figure1,...
    'Position',[0.55 0.22 0.33 0.3]);
hold(axes4,'on');

scatter(subm_all(~ind & T_inv),vcd_aod_all(~ind & T_inv),50,'k','filled','marker','pentagram');
scatter(subm_all(ind & T_inv),vcd_aod_all(ind & T_inv),50,wspd_all(ind & T_inv),'filled','marker','pentagram');

scatter(subm_all(~ind & ~T_inv),vcd_aod_all(~ind & ~T_inv),25,'k','filled');
scatter(subm_all(ind & ~T_inv),vcd_aod_all(ind & ~T_inv),25,wspd_all(ind & ~T_inv),'filled');

ylim([0,12e13])
% xlim([2,12]*1000)

box on
% Create xlabel
xlabel('0.1\mum < D_p < 0.5\mum (cm^-^3)');

% Create colorbar
c=colorbar('peer',axes4,'location','southoutside','Position',...
    [0.3333 0.06 0.33333 0.03]);
ylabel(c,'Windspeed (m/s)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % %% plot separate figures by aerosol size
% % 
% % figure1 = figure;
% % colormap(jet(14));
% % 
% % %% plot upper half of VCD vs supermicrom particle conc.
% % axes1 = axes('Parent',figure1,...
% %     'Position',[0.13 0.3 0.33 0.6]);
% % hold(axes1,'on');
% % 
% % % scatter(sm_all(~ind),vcd_top_all(~ind),25,'k','filled');
% % % scatter(sm_all(ind),vcd_top_all(ind),25,wspd_all(ind),'filled');
% % 
% % % with temperature inversion below PEARL
% % scatter(sm_all(~ind & T_inv),vcd_top_all(~ind & T_inv),50,'k','filled','marker','pentagram');
% % scatter(sm_all(ind & T_inv),vcd_top_all(ind & T_inv),50,wspd_all(ind & T_inv),'filled','marker','pentagram');
% % 
% % % no temperature inversion below PEARL
% % scatter(sm_all(~ind & ~T_inv),vcd_top_all(~ind & ~T_inv),25,'k','filled');
% % scatter(sm_all(ind & ~T_inv),vcd_top_all(ind & ~T_inv),25,wspd_all(ind & ~T_inv),'filled');
% % 
% % ylim([0,4.5e13])
% % 
% % box on
% % % Create ylabel
% % ylabel('BrO VCD_{0.7-4km} (mol/cm^2)');
% % xlabel('D_p > 1\mum (cm^-^3)');
% % 
% % %% plot full VCD vs supermicrom particle conc.
% % axes2 = axes('Parent',figure1,...
% %     'Position',[0.55 0.3 0.33 0.6]);
% % hold(axes2,'on');
% % 
% % scatter(sm_all(~ind & T_inv),vcd_aod_all(~ind & T_inv),50,'k','filled','marker','pentagram');
% % scatter(sm_all(ind & T_inv),vcd_aod_all(ind & T_inv),50,wspd_all(ind & T_inv),'filled','marker','pentagram');
% % 
% % scatter(sm_all(~ind & ~T_inv),vcd_aod_all(~ind & ~T_inv),25,'k','filled');
% % scatter(sm_all(ind & ~T_inv),vcd_aod_all(ind & ~T_inv),25,wspd_all(ind & ~T_inv),'filled');
% % 
% % ylim([0,12e13])
% % 
% % box on
% % % Create xlabel
% % ylabel('BrO VCD_{0-4km} (mol/cm^2)');
% % xlabel('D_p > 1\mum (cm^-^3)');
% % 
% % % Create colorbar
% % c=colorbar('peer',axes2,'location','southoutside','Position',...
% %     [0.3333 0.14 0.33333 0.03]);
% % ylabel(c,'Windspeed (m/s)')
% % 



% % % figure1 = figure;
% % % colormap(jet(14));
% % % 
% % % 
% % % %% plot upper half of VCD vs submicron particle conc.
% % % axes1 = axes('Parent',figure1,...
% % %     'Position',[0.13 0.3 0.33 0.6]);
% % % hold(axes1,'on');
% % % 
% % % scatter(subm_all(~ind & T_inv),vcd_top_all(~ind & T_inv),50,'k','filled','marker','pentagram');
% % % scatter(subm_all(ind & T_inv),vcd_top_all(ind & T_inv),50,wspd_all(ind & T_inv),'filled','marker','pentagram');
% % % 
% % % scatter(subm_all(~ind & ~T_inv),vcd_top_all(~ind & ~T_inv),25,'k','filled');
% % % scatter(subm_all(ind & ~T_inv),vcd_top_all(ind & ~T_inv),25,wspd_all(ind & ~T_inv),'filled');
% % % 
% % % ylim([0,4.5e13])
% % % % xlim([2,12]*1000)
% % % 
% % % xlabel('0.1\mum < D_p < 0.5\mum (cm^-^3)');
% % % ylabel('BrO VCD_{0.7-4km} (mol/cm^2)');
% % % 
% % % box on
% % % 
% % % %% plot full VCD submicron particle conc.
% % % axes2 = axes('Parent',figure1,...
% % %     'Position',[0.55 0.3 0.33 0.6]);
% % % hold(axes2,'on');
% % % 
% % % scatter(subm_all(~ind & T_inv),vcd_aod_all(~ind & T_inv),50,'k','filled','marker','pentagram');
% % % scatter(subm_all(ind & T_inv),vcd_aod_all(ind & T_inv),50,wspd_all(ind & T_inv),'filled','marker','pentagram');
% % % 
% % % scatter(subm_all(~ind & ~T_inv),vcd_aod_all(~ind & ~T_inv),25,'k','filled');
% % % scatter(subm_all(ind & ~T_inv),vcd_aod_all(ind & ~T_inv),25,wspd_all(ind & ~T_inv),'filled');
% % % 
% % % ylim([0,12e13])
% % % % xlim([2,12]*1000)
% % % 
% % % box on
% % % % Create xlabel
% % % xlabel('0.1\mum < D_p < 0.5\mum (cm^-^3)');
% % % ylabel('BrO VCD_{0-4km} (mol/cm^2)');
% % % 
% % % % Create colorbar
% % % c=colorbar('peer',axes2,'location','southoutside','Position',...
% % %     [0.3333 0.14 0.33333 0.03]);
% % % ylabel(c,'Windspeed (m/s)')
