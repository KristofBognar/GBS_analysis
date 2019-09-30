% get cross-sections for peaks found in low signal spectra
% wavelengths are obtained by doing a calibration on the nearest proper
% spectrum -- the standard cal file might not be accurate due to wavelength
% shifts at night


load('/home/kristof/work/documents/GBS_reanalysis_update_sep17/PEARL-GBS/saved.mat');


% plot(p50(:,1),p50(:,2)), hold on
% plot(cal273,p273(:,2))
% plot(p296(:,1),p296(:,2))


% plot(cal286,p2861(:,2)), hold on
% plot(cal286,p2862(:,2))
% plot(p53(:,1),p53(:,2)*5), hold on
% plot(p54(:,1),p54(:,2))


%% peak at 436 nm
% use d273 for double peak at 550nm and peak at 436nm

% since sza is not high enough, solar spectrum contaminates the lamp peaks
% scale kurucz spectrum to each peak manualy to subtract the background
% plot(cal273,p273(:,2)), hold on 
% plot(cal273,kur273*mean(p273(:,2)./kur273)*0.8)
% xlim([432,440])

spec=p273(:,2)-kur273*mean(p273(:,2)./kur273)*0.8;

plim=[434.5,437];

pind=find(cal273>=plim(1) & cal273<=plim(2));

peak=spec(pind);
lambda=cal273(pind);

% out=[lambda,peak];
% dlmwrite(['X436nm_P1_2007_273.xs'],out,'\t');

%% double peak (~545 nm)

% plot(cal273,p273(:,2)), hold on 
% plot(cal273,kur273*mean(p273(:,2)./kur273)*1.25)
% xlim([500,555])

spec=p273(:,2)-kur273*mean(p273(:,2)./kur273)*1.25;

plim=[530,560];

pind=find(cal273>=plim(1) & cal273<=plim(2));

peak=spec(pind);
lambda=cal273(pind);

peak(194:195)=[];
lambda(194:195)=[];

% plot(lambda,peak), hold on
% peak=boxcar(lambda,peak,2)';
% plot(lambda,peak)
% out=[lambda,peak];
% dlmwrite(['X545nm_P1_2007_273.xs'],out,'\t');


%% peak at 498 nm
% use d273 for double peak at 550nm and peak at 436nm

% plot(cal273,p273(:,2)), hold on 
% plot(cal273,kur273*mean(p273(:,2)./kur273)*1.02)
% xlim([470,526])

spec=p273(:,2)-kur273*mean(p273(:,2)./kur273)*1.02;

plim=[494,504];

pind=find(cal273>=plim(1) & cal273<=plim(2));

peak=spec(pind);
lambda=cal273(pind);

% plot(lambda,peak), hold on
% peak=boxcar(lambda,peak,2)';
% plot(lambda,peak)
% out=[lambda,peak];
% dlmwrite(['X498nm_P1_2007_273.xs'],out,'\t');

%% broad peak at ~465 nm
% doesn't seem appear in the residuals

% plot(cal286,p2861(:,2)), hold on
% plot(cal286,p2862(:,2))
m=(p2861(:,2)+p2862(:,2))./2;
plot(m)
ind_out=[773,774,806,832,833,836:838,913,914,1023:1027,1209,1337:1343,...
         1371,1657:1660,1664:1667,16891695:1697];

% remove peak at 498 nm, it's fitted separately
m(1511:1541)=m(1511:1541)-5*exp(-((cal286(1511:1541)-498)./1.942).^2);
     
m(ind_out)=[];
lambda_tmp=cal286;
lambda_tmp(ind_out)=[];

% plot(lambda_tmp,m)
% xlim([430,510])


plim=[429,510];

pind=find(lambda_tmp>=plim(1) & lambda_tmp<=plim(2));

peak=m(pind);
lambda=lambda_tmp(pind);

plot(lambda,peak), hold on
peak=boxcar(lambda,peak,15)';
plot(lambda,peak)
% out=[lambda,peak];
% dlmwrite(['X465nm_P1_2007_286.xs'],out,'\t');



