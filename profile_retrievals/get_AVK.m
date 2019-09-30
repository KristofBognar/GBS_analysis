% read averaging kernel data from current folder

% make list of all files
temp = dir('*.dat'); 
f_list = {temp.name}; % cell array of file names

%% plot single AVK

avk=dlmread(f_list{12},'',1,5);
alt=0:0.2:4;

figure(1)

for i=1:20
    
    plot(avk(:,i),alt,'linewidth',1.3), hold on
    
end

legend('0.0 km','0.2 km','0.4 km','0.6 km','0.8 km','1.0 km','1.2 km','1.4 km',...
       '1.6 km','1.8 km','2.0 km','2.2 km','2.4 km','2.6 km','2.8 km',...
       '3.0 km','3.2 km','3.4 km','3.6 km','3.8 km','location','eastoutside')

xlabel('Averaging kernel')
ylabel('Altitude (km)')

%% save diagonal of each avk matrix

% diag_all=[];
% for i=1:length(f_list)
%     
%     tmp=dlmread(f_list{i},'',1,5);
%     diag_all=[diag_all, diag(tmp)];
%     sum(tmp,1)
% end
% bottom=sum(diag_all(1:3,:),1);
% top=sum(diag_all(5:end,:),1);
% layers=[bottom;diag_all(4,:);top];

