% function [ output_args ] = HOBr_reaction_rate( input_args )
%
%
% debug
year=2017;


%% load OPC data

if year==2016
    load(['/home/kristof/work/SMPS/opc_size_dist_0316.mat']);
elseif year==2017
    load(['/home/kristof/work/SMPS/opc_size_dist_0317.mat']);
    load(['/home/kristof/work/weather_stations/ridge_lab/PWS_03_2017.mat']);
    
end

% particle size bins from OPC file, in microns (for last bin only the lower
% bound is known, but no particles are that large anyway)
sizes=[0.4;0.75;1.5;3.5;7.5;15];

% bins to include
b1=3;
b2=4;

% convert times to fractional time
[ft_opc]=fracdate(time,'dd-eee-yyyy hh:mm:ss');

% get pressure, temperature at ridge lab altitude
T_smooth=boxcar(ft_wnd,T+273.15,5)';
P_smooth=boxcar(ft_wnd,P,5)';

T_RL=interp1(ft_wnd,T_smooth,ft_opc);
P_RL=interp1(ft_wnd,P_smooth,ft_opc);

% get air number density for ppt calculation (cm-3)
num_dens=((6.022141e23*P_RL)./(8.3144598*T_RL))*1e-6; 


ft_opc=ft_opc-58; % convert to days of march, 2017

% convert dN/dlog10(Dp) to N
Dp=[Dp,20000];
logDp=log10(Dp(2:end)./Dp(1:end-1));
for i=1:length(logDp)
    Dp_data(:,i)=Dp_data(:,i)*logDp(i);
end

% get supermicron particle info (cm-3)
supermicron_full=sum(Dp_data(:,b1:b2),2);

% get mean diameter in microns
d_mean=(Dp_data(:,b1:b2)*sizes(b1:b2))./supermicron_full;

% filter out bad data
ind=find(isnan(T_RL) | supermicron_full>400 | supermicron_full==0);
ft_opc(ind)=[];
T_RL(ind)=[];
P_RL(ind)=[];
num_dens(ind)=[];
Dp_data(ind,:)=[];


% mean molecular speed in cm/s (!!)
m=0.097/6.022e23; % mass of HOBr molecule in kg
kb=1.38065e-23; % boltzmann constant in J/K
% T=250; % temperature in K, SHOULD USE REAL DATA!!
c_mean=sqrt(2 * (3*kb*T_RL/2) / m) *100; % speed, converted to cm/s

K_tot=zeros(length(ft_opc),6);
br_tot=zeros(length(ft_opc),6);

for i=b1:b2
    %% get surface area concentration (cm2/cm3) and volume (cm3/cm3)
    
    Area=4*pi*(sizes(i)*1e-4).^2 * Dp_data(:,i);
    Volume=(4/3)*pi*(sizes(i)*1e-4).^3 * Dp_data(:,i);
    
    %% Estimate heterogeneous HOBr reaction rate (Tang et al., 2014)

    % Knudsen number for diffusion limitation
    % P=1e5; % pressure in Pa, SHOULD USE REAL DATA!!
    diff_coeff=84./(P_RL*0.00750062); % converted from torr cm2 s-1 to cm2 s-1, value of 84 from Tang et al
    Kn=(6*diff_coeff)./(c_mean * sizes(i)*1e-4);

    % correction for diffusion limitation for large particles
    correction=(0.75+0.286*Kn)./(Kn.*(Kn+1));

    % uptake coefficient (Roberts et al., 2014)
    gamma=0.6; % gamma value for small particles
    gamma_eff=(1/gamma + correction).^-1;

    % reaction rate in s-1
    K=gamma_eff.*c_mean.*(Area)./4;
    
    % add to total reaction rate
    K_tot(:,i)=K;


    %% monthly mean Br conc. in aerosol (Alert, 1980-86)
    % comparable to observed BrO concentrations
    % br_mean=30*1e-9 * 1e-6; % nanograms/m3, converted to g/cm3
    % br_mean=br_mean*6.022e23/80; % converted to n.o. molecules/cm3

    %% Br- concentration (or mass) from seawater ratios (Quinn et al., 2002)
    % cl-/br- ratio -- varies a lot in aerosol/snow (Simpson et al., 2007)
    % cl_to_br=590.91; 

%     % alternatively, specify Br- enhancement (in %) over SW ratio
%     enhance=10;
%     cl_to_br=1/((enhance/65000)+(1/650));
% 
%     mass_tot=Volume*2.17; % total aerosol mass in g/cm3 (g per cm3 of air, NOT density)
%     br_mass_tot=mass_tot/(650*3.65); % br- mass using seawater cl/br ratio
%     br_mass_free=mass_tot/(cl_to_br*3.65)-br_mass_tot; % mass on top of seawater ratio 
%     br_conc=br_mass_free*6.022e23/80; % concentration in molecules/cm3

    % fix percent of sea water ratio released
    mass_tot=Volume*2.17; % total aerosol mass in g/cm3 (g per cm3 of air, NOT density)
    br_mass_tot=(mass_tot/(650*3.65))/2; % br- mass using seawater cl/br ratio, and assuming 50% SSA content
    br_conc=(br_mass_tot*6.022e23/80)*0.03; % concentration in molecules/cm3, with only some fraction released
    
    % add to total concentration
    br_tot(:,i)=br_conc;
    

    % % % load(['/home/kristof/work/SMPS/size_dist_0317.mat'])
    % % % 
    % % % Dp=[Dp,500];
    % % % logDp=log10(Dp(2:end))-log10(Dp(1:end-1));
    % % % 
    % % % N=Dp_data*logDp';


end



figure(1)

for i=b1:b2
    
    figure(1)
    plot(ft_opc,K_tot(:,i)), hold on
    
    figure(2)
    plot(ft_opc,br_tot(:,i)*1e12./num_dens), hold on
    
end

figure(1)

plot([ft_opc(1),ft_opc(end)],[0.0025,0.0025]), hold on % Peterson et al, 2017 supplement

xlabel(['Fractional date ' num2str(year) ' (UTC)'])
ylabel('HOBr het. reaction rate (s^-^1)')

xlim([7,21])

legend('D_p=1.5 \mum','D_p=3.5 \mum','Approx. HOBr photolysis rate')

figure(2)

legend('D_p=1.5 \mum','D_p=3.5 \mum')

xlabel(['Fractional date ' num2str(year) ' (UTC)'])
% ylabel('[Br^-] (molec/cm3) on top of SS Cl^-/Br^- ratio')
ylabel('[Br^-] (pptv) on top of SS Cl^-/Br^- ratio')
% title(['Cl^-/Br^- = ' num2str(cl_to_br)])
xlim([7,21])

figure(99)
hold on
h99_2=plot(ft_opc,sum(br_tot(:,b1:b2),2)*1e12./num_dens); hold on
xlabel('Days of March, 2017 (UTC)')
ylabel('BrO at 600 m (pptv)')

% legend([h99_2,h99_1],{'Estimate','MAX-DOAS'},'Location','NorthWest')


xlim([7,21])

