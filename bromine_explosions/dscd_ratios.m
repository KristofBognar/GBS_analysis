function [ ratio_all, time_all ] = dscd_ratios( ea_in, dscd_only )
% [ ratio_all, time_all ] = dscd_ratios( ea_in )
%   Calculates O4-normalised ratio of -1deg elevation angle to another
%   angle, specified in ea_in
%       Ratio is calculated only if both measurements are in the same scan
%   Also outputs the mean time for the ratio (mean of measurememt times at each elev angle)
%
%   if dscd_only=1, the code returns the dscd for the specified elevation
%   angle (normalized by O4)
%
%@Kristof Bognar, June 2020

% input check
if ~any([0,1]==ea_in), error('add support for other elevation angles'), end

if nargin==1, dscd_only=0; end

%% load data

data=[];
load('/home/kristof/aurora/ground/eureka/gbs/pearl-gbs/2016/output/maxdoas_d66-131.mat')
data=maxdoas;
load('/home/kristof/aurora/ground/eureka/gbs/pearl-gbs/2017/QDOAS_output/maxdoas_d66-94.mat')
data=[data;maxdoas];
load('/home/kristof/aurora/ground/eureka/gbs/pearl-gbs/2018/QDOAS_output/maxdoas_d64-151.mat')
data=[data;maxdoas];
load('/home/kristof/aurora/ground/eureka/gbs/pearl-gbs/2019/QDOAS_output/maxdoas_d64-151.mat')
data=[data;maxdoas];

%% get rid of missing times
% matlab import occasionally converts valid datatime to NaT for some reason
ind_nat=isnat(data.DateTime);
ind_ok=~isnat(data.DateTime);

data.DateTime(ind_nat)=...
         ft_to_date(data.Fractionalday(ind_nat)-1,data.Year(ind_nat));

if sum(isnat(data.DateTime)), error('Could not fill all NaTs'); end

%% get output

if dscd_only==0
    %% get dSCD ratio

    % change -1deg to create unique differences (1-2, 0-1, -1-0 would all be -1)
    data.Elevviewingangle(data.Elevviewingangle==-1)=101; % arbitrary number

    % remove 0deg lines if interested in 1deg
    if ea_in==1, data(data.Elevviewingangle==0,:)=[]; end

    % differences where given elev angle is followed by -1deg are easy to identify
    diffs=[0;diff(data.Elevviewingangle)];

    % indices of -1deg, only when preceeded by given ea
    if ea_in==1
        indm1=find(diffs==100);
    elseif ea_in==0
        indm1=find(diffs==101);
    end

    % ideices of given ea
    ind_ea=indm1-1;

    % scale by O4
    dscd_corr=data.BrOSlColbro./data.O4SlColo4;

    % get ratio and mean time
    ratio_all=dscd_corr(indm1)./dscd_corr(ind_ea);
    time_all=mean([data.DateTime(indm1),data.DateTime(ind_ea)],2);

else
    %% dscd value only
    
    ratio_all=data.BrOSlColbro(data.Elevviewingangle==ea_in)./...
              data.O4SlColo4(data.Elevviewingangle==ea_in);
    time_all=data.DateTime(data.Elevviewingangle==ea_in);
    
end

end



%% leftover bits of code for testing/plotting
% % elevation angles in the file (might change from year to year)
% elevs=unique(maxdoas.Elevviewingangle);
% elevs=setdiff(elevs,0);
% 
% % colors for plotting
% colors={'r.','g.','b.','y.','m.','c.','k.','ko'};
% 
% msize=16;
% 
% % empty legend cell
% legends=cell(1,length(elevs));
% 
% % initialize figure
% figure(99)
% set(gcf, 'Position', [100, 100, 900, 600]);
% 
% for i=1:length(elevs)
% 
%     ax_bro=subplot(211);
%     
%     % find indices with current elevation
%     ind=(data.Elevviewingangle==elevs(i) & data.BrORMS<0.003);
%          
%     % plot
%     if i==8, msize_tmp=4; else msize_tmp=msize; end
%     plot(data.DateTime(ind),data.BrOSlColbro(ind),colors{i},'markersize', msize_tmp); hold on
%         
%     % create legend entry
%     legends{i}=[num2str(elevs(i)) '\circ'];
%             
% end 
% 
% legend(legends,'location','northeast','Orientation','horizontal')
% 
% for i=1:length(elevs)
% 
%     ax_o4=subplot(212);
%     
%     % find indices with current elevation
%     ind=(data.Elevviewingangle==elevs(i) & data.O4RMS<0.003);
%          
%     % plot
%     if i==8, msize_tmp=4; else msize_tmp=msize; end
%     plot(data.DateTime(ind),data.O4SlColo4(ind),colors{i},'markersize', msize_tmp); hold on
%         
%     % create legend entry
%     legends{i}=[num2str(elevs(i)) '\circ'];
%             
% end 
% 
% linkaxes([ax_bro,ax_o4],'x');

