function get_sonde_for_profile(start_yy,start_mm,start_dd,end_yy,end_mm,end_dd)
% create interpolated radiosonde profiles for middle of measurement day to
% be used in profile retrieval
%
% INPUT:    year, month and day of first and last day sonde data is required for
%
% OUTPUT:   Formatted radiosonde and ozonesonde files in specified directory

%% input parameters
% time set to local noon
t1 = datetime(start_yy,start_mm,start_dd,18,0,0);
t2 = datetime(end_yy,end_mm,end_dd,18,0,0);
dates = t1:t2;

% calculate fractional time of measurements
[ft,year]=fracdate(dates);

if length(year)~=length(ft)
    year=ones(size(ft))*year;
end
    

savedir='/home/kristof/work/profile_retrievals/data_for_profile_retrieval/sonde_data/';


%% Loop over times

for ii=1:length(dates)
    %% radiosonde data

    % interpolate radiosonde data to given time
    [z,Pa,K]=interp_radiosonde(year(ii),ft(ii));
    asl_km=z/1000;
    hPa=Pa/100;

    % interpolate to 100m grid
    grid=[0:0.1:50];
    hPa=interp1(asl_km,hPa,grid)';
    hPa(1)=Pa(1)/100;
    k1=K(1);
    K=interp1(asl_km,K,grid)';
    K(1)=k1;
    asl_km=grid';

    % create table from data
    radiosonde=table(asl_km,hPa,K);

    % remove NaNs
    ind=find(isnan(hPa));

    if isempty(ind) % do nothing
    elseif ind==find(isnan(K))
        radiosonde=radiosonde(1:ind(1)-1,:);
    else
        error('P and T profiles end at different altitudes') % should not be possible
    end

    % create filename and write file
    savename=[savedir, 'radiosonde_', datestr(dates(ii),'yymmdd'), '.dat'];
    writetable(radiosonde,savename,'Delimiter',',');

    %% ozonesonde data

    % interpolate ozonesonde data to given time
    [z,vmr,dens,~]=interp_ozonesonde(year(ii),ft(ii));
    alt=z/1000;
    o3_conc=vmr.*dens;

    % interpolate to 100m grid
    k1=o3_conc(1);
    o3_conc=interp1(alt,o3_conc,grid)';
    o3_conc(1)=k1;
    alt=grid';

    % create table from data
    ozonesonde=table(alt,o3_conc);

    % remove NaNs
    ind=find(isnan(o3_conc));
    ozonesonde=ozonesonde(1:ind(1)-1,:);

    % create filename and write file
    savename=[savedir, 'ozonesonde_', datestr(dates(ii),'yymmdd'), '.dat'];
    writetable(ozonesonde,savename,'Delimiter',',','WriteVariableNames',false);

end

end
