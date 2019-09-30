function [ ] = plot_opc( )
%PLOT_OPC Summary of this function goes here
%   Detailed explanation goes here

year=2017;

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

% filter out bad data
ind=find(supermicron_full==0 | supermicron_full>400);
supermicron_full(ind)=[];
ft_opc(ind)=[];

if year==2017
    ind=[7542;7543];
    supermicron_full(ind)=[];
    ft_opc(ind)=[];
end

plot(ft_opc-58,supermicron_full,'b-')
grid on

xlim([7,22.1])

ylabel('Dp > 1\mum (cm^-^3)')
hold on

end

