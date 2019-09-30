function [ time_out ] = time_from_SZA( year, doy, sza, loc )
%TIME_FROM_SZA(year, doy, sza, loc) 
%   Function to calculate approximate time from date, location and solar
%   zenith angle data
%   Time resolution is hard-coded in the function, see variable 'step'
%
%   INPUT:  year: either single value, or array with the same size as SZA and doy
%           doy: day of the year, jan. 1 = 1, same size as sza
%           sza: array of solar zenith angle values
%           loc: location identifier (e.g. 'E' for Eureka)
%
%   OUTPUT: ref_time_out: the approximate time (matlab serial time, UTC) of each SZA value. Since
%               each SZA repeats twide a day, two values are reported.
%               first column: closest time
%               second column: second closest time, on the opposite side of
%                  noon by default
%               If the closest time is +- 30 min from 18:00 (approximate
%               local noon) then only the first value is reported

%% Location data
if strcmp(loc,'E')
    lon = -86.416;
    lat = 80.053;
    alt = 0.61;
else
    error('Location not recognized')
end

%% get date (in datetime format) from doy
date=datetime(yeartime(year)+doy-1,'convertfrom','datenum');

% convert inputs to column vector if necessary
if size(date,2)~=1, date=date'; end
if size(sza,2)~=1, sza=sza'; end

%% get solar elevation values for given dates in given intervals

% time step for getting SZA values (determines how close the time is to actual SZA)
% in hours
step=15/60;

% time values to loop through
hours=[0:step:23];

% array to store solar elevations
el=NaN(length(doy),length(hours));

% get solar elevations for each hour
for i=1:length(hours)
    [~,el(:,i)]=SolarAzEl(date+hours(i)/24,lat,lon,alt);
end

%% find closest times to SZA

% repeat SZA vector, one column for each time, and convert to elevation
input_el=90-repmat(sza,1,length(hours));

% sort the differences between hourly elev and input SZA
[~,sortind]=sort(abs(el-input_el),2);

% times of SZA are the two minima from el-input_el (first two indices in sortind)
time=NaN(length(sza),2);
for i=1:length(sza), time(i,:)=hours(sortind(i,1:2)); end

% chatch missing input
ind=find(sza==0);
time(ind,:)=NaN(length(ind),2);

%% format output

% find times close to noon
ind=find(abs(time(:,1)-18)<=.5);
time(ind,2)=NaN(length(ind),1);

% % convert to date number
time_out=datenum(date+time(:,1)/24);
time_out=[time_out, datenum(date+time(:,2)/24)];

% time_out=time;

end

