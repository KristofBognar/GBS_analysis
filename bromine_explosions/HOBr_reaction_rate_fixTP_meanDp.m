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
end

% particle size bins from OPC file, in microns (for last bin only the lower
% bound is known, but no particles are that large anyway)
sizes=[0.4;0.75;1.5;3.5;7.5;15];

% bins to include
b1=3;
b2=5;

% convert times to fractional time
[ft_opc]=fracdate(time,'dd-eee-yyyy hh:mm:ss');

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

% get total surface area concentration (cm2/cm3) and volume (cm3/cm3)
A_tot=0;
V_tot=0;
for i=b1:b2
    A_tot=A_tot + 4*pi*(sizes(i)*1e-4).^2 * Dp_data(:,i);
    V_tot=V_tot + (4/3)*pi*(sizes(i)*1e-4).^3 * Dp_data(:,i);
end

% filter out bad data
ind=find(supermicron_full==0 | supermicron_full>400);
supermicron_full(ind)=[];
ft_opc(ind)=[];
d_mean(ind)=[];
A_tot(ind)=[];
V_tot(ind)=[];


%% Estimate heterogeneous HOBr reaction rate (Tang et al., 2014)

% mean molecular speed in cm/s (!!)
m=0.097/6.022e23; % mass of HOBr molecule in kg
kb=1.38065e-23; % boltzmann constant in J/K
T=250; % temperature in K, SHOULD USE REAL DATA!!
c_mean=sqrt(2 * (3*kb*T/2) / m) *100; % speed, converted to cm/s


% Knudsen number for diffusion limitation
P=1e5; % pressure in Pa, SHOULD USE REAL DATA!!
diff_coeff=84/(P*0.00750062); % converted from torr cm2 s-1 to cm2 s-1, value of 84 from Tang et al
Kn=(6*diff_coeff)./(c_mean * d_mean*1e-4);

% correction for diffusion limitation for large particles
correction=(0.75+0.286*Kn)./(Kn.*(Kn+1));

% uptake coefficient
gamma=0.6; % gamma value for small particles
gamma_eff=(1/gamma + correction).^-1;

% reaction rate in s-1
K=gamma_eff.*c_mean.*(A_tot)./4;

figure()
plot([ft_opc(1),ft_opc(end)],[0.0025,0.0025]), hold on % Peterson et al, 2017 supplement
plot(ft_opc,K), hold on

xlabel(['Fractional date ' num2str(year) ' (UTC)'])
ylabel('HOBr het. reaction rate (s^-^1)')
legend('Approx. HOBr photolysis rate')
xlim([7,21])

%% monthly mean Br conc. in aerosol (Alert, 1980-86)
% comparable to observed BrO concentrations
% br_mean=30*1e-9 * 1e-6; % nanograms/m3, converted to g/cm3
% br_mean=br_mean*6.022e23/80; % converted to n.o. molecules/cm3

%% Br- concentration (or mass) from seawater ratios (Quinn et al., 2002)
% cl-/br- ratio -- varies a lot in aerosol/snow (Simpson et al., 2007)
% cl_to_br=590.91; 

% alternatively, specify Br- enhancement (in %) over SW ratio
enhance=10;
cl_to_br=1/((enhance/65000)+(1/650));

mass_tot=V_tot*2.17; % total aerosol mass in g/cm3 (g per cm3 of air, NOT density)
br_mass_tot=mass_tot/(650*3.65); % br- mass using seawater cl/br ratio
br_mass_free=mass_tot/(cl_to_br*3.65)-br_mass_tot; % mass on top of seawater ratio 
br_conc=br_mass_free*6.022e23/80; % concentration in molecules/cm3

figure()
plot(ft_opc,br_conc)

xlabel(['Fractional date ' num2str(year) ' (UTC)'])
ylabel('[Br^-] (molec/cm3) on top of SS Cl^-/Br^- ratio')
title(['Cl^-/Br^- = ' num2str(cl_to_br)])
xlim([7,21])

% % % load(['/home/kristof/work/SMPS/size_dist_0317.mat'])
% % % 
% % % Dp=[Dp,500];
% % % logDp=log10(Dp(2:end))-log10(Dp(1:end-1));
% % % 
% % % N=Dp_data*logDp';


% end

