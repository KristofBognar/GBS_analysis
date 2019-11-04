% plot VCDs

figure(1)

subplot(211)

% VCD plot
gscatter(VCD_table.fd,VCD_table.mean_vcd,VCD_table.ampm,'br','..'); hold on,
% legend('UT-GBS a.m.','UT-GBS p.m.');
gscatter(VCD_table_P_no2.fd,VCD_table_P_no2.mean_vcd,VCD_table_P_no2.ampm,'br','oo');
legend('UT-GBS a.m.','UT-GBS p.m.','PEARL-GBS a.m.','PEARL-GBS p.m.','location','northwest');
ylabel('NO_2 molecule/cm^2')
plot([83.66,83.66],[0,15]*10^15, 'k-')
plot([183.66,183.66],[0,15]*10^15, 'k-')

% title('VCDs')
xlim([40,274])
grid on

subplot(212)
i=1;
maxind=size(VCD_table.day,1);

% Calculate daily pm-am difference if there are 2 measurements for the day
diff_u=[];
day_u=[];
for i=49:274
    k=find(VCD_table.day(:,1)==i);
    if size(k,1) == 2
        diff_u=[diff_u, VCD_table.mean_vcd(k(2),1) - VCD_table.mean_vcd(k(1),1)];
        day_u=[day_u, VCD_table.day(k(1),1)];
    end
end

% % calcluate mean by averaging n days, going backwards
% n=5;
mday_u=[];
mdiff_u=[];
% for i=size(day_u,2):-n:mod(size(day_u,2),n)+n
%     mday_u=[mday_u, mean(day_u(i-n+1:i))];
%     mdiff_u=[mdiff_u, mean(diff_u(i-n+1:i))];
% end
% mday_u=[mday_u,mean(day_u(1:i-n))];
% mdiff_u=[mdiff_u,mean(diff_u(1:i-n))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same for PGBS data
i=1;
maxind=size(VCD_table_P_no2.day,1);

% Calculate daily pm-am difference
diff_p=[];
day_p=[];
for i=50:129
    k=find(VCD_table_P_no2.day(:,1)==i);
    if size(k,1) == 2
        diff_p=[diff_p, VCD_table_P_no2.mean_vcd(k(2),1) - VCD_table_P_no2.mean_vcd(k(1),1)];
        day_p=[day_p, VCD_table_P_no2.day(k(1),1)];
    end
end

% % calcluate mean by averaging n days, going backwards
mday_p=[];
mdiff_p=[];
% for i=size(day_p,2):-n:mod(size(day_p,2),n)+n
%     mday_p=[mday_p, mean(day_p(i-n+1:i))];
%     mdiff_p=[mdiff_p, mean(diff_p(i-n+1:i))];
% end
% mday_p=[mday_p,mean(day_p(1:i-n))];
% mdiff_p=[mdiff_p,mean(diff_p(1:i-n))];

% calculate 5-day mean here for both instruments
for i=40:5:275
    
    % UT-GBS
    % find indices for given day, making sure there are numbers to average
    ind1=find(day_u==(i-2));
    temp=1;
    % move in from boundary if it's missing
    while isempty(ind1) && temp < 3
        ind1=find(day_u==(i-2+temp)); 
        temp=temp+1;        
    end 
 
    % same with other index
    ind2=find(day_u==(i+2));
    temp=1;
    while isempty(ind2) && temp < 3
        ind2=find(day_u==(i+2-temp)); 
        temp=temp+1;        
    end 

    % make sure there's data if there's only one datapoint in bin
    if isempty(ind1), ind1=ind2; end
    if isempty(ind2), ind2=ind1; end
    
    nextval=mean(diff_u(ind1:ind2));
    mdiff_u=[mdiff_u, nextval];
     
    if ~isempty(nextval),mday_u=[mday_u, i];end
    
    % PEARL-GBS
    % find indices for given day
    ind1=find(day_p==(i-2));
    temp=1;
    while isempty(ind1) && temp < 3
        ind1=find(day_p==(i-2+temp)); 
        temp=temp+1;        
    end 
 
    ind2=find(day_p==(i+2));
    temp=1;
    while isempty(ind2) && temp < 3
        ind2=find(day_p==(i+2-temp)); 
        temp=temp+1;        
    end 

    if isempty(ind1), ind1=ind2; end
    if isempty(ind2), ind2=ind1; end
    
    nextval=mean(diff_p(ind1:ind2));
    mdiff_p=[mdiff_p, nextval];
    if ~isempty(nextval), mday_p=[mday_p, i];end
    
end


plot(day_u, diff_u, 'r.','markersize',14), hold on,
plot(day_p, diff_p, 'b.','markersize',14), hold on,

plot(mday_u, mdiff_u,'r-', 'linewidth',2), hold on
plot(mday_p, mdiff_p,'b-', 'linewidth',2), hold on

legend('UT-GBS','PEARL-GBS','UT-GBS avg','PEARL-GBS avg','location','best')
ylim([-5e14,10e14])

plot([83.66,83.66],[-5,10]*10^14, 'k-')
plot([183.66,183.66],[-5,10]*10^14, 'k-')


xlabel('Day of the year, 2016 (UTC)')
ylabel('\Delta NO_2 molecule/cm^2')
% title('Diurnal variation')
xlim([40,274])
grid on


