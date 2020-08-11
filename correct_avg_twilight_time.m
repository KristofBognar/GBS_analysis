function [ table_out ] = correct_avg_twilight_time( table_in )
%[ table_out ] = correct_avg_twilight_time( table_in )
%
% Function to fix bug in Preprocs_ASC.m for averaged twilight spectra:
%
%       For SZA bins that include measurements around UTC midnight, the
%       time average was calculated as arithmetic mean, not proper datatime
%       average. As a result, at most one averaged spectrum per L2 file has
%       the wrong time information.
%       This doesn't affect processing in QDOAS, as the spectrum is
%       included in the correct daily file (i.e. it's processed with the
%       correct daily reference).
%
% The code replaces the bad values in the following fields:
%   Fractionalday, Fractionaltime, SZA, SolarAzimithAngle, DateDDMMYYYY,
%   Timehhmmss
%
% Returns the full table, with the bad rows fixed
%
%@Kristof Bognar, Aug. 2020

table_out=table_in;

%% find bad values
% all bugged spectra have exeriment time<0, due to the day change not being
% considered in the time averaging
bad_ind=find(table_out.TotalExperimentTimesec<0);

% if no bad values, stop without changing anythin
if isempty(bad_ind), return, end

% add date to time field, for correct averages and conversion to fractional time
table_out.Timehhmmss=table_out.DateDDMMYYYY+timeofday(table_out.Timehhmmss);
table_out.Timehhmmss.Format='HH:mm:ss';

% replace all bad dates/times with NaT
table_out.DateDDMMYYYY(bad_ind)=NaT;
table_out.Timehhmmss(bad_ind)=NaT;


%% fill bad values
for i=bad_ind'

    % check if average or extrapolation is needed
    if table_out.SpecNo(i+1)-table_out.SpecNo(i)~=1
        
        % bad line is the last spectrum of the day; extrapolate from
        % previous two measurements
        table_out.Timehhmmss(i)=table_out.Timehhmmss(i-1)+diff(table_out.Timehhmmss(i-2:i-1));
        
        % SZA, SAA also needs extrapolation
        table_out.SolarAzimuthAngle(i)=...
            table_out.SolarAzimuthAngle(i-1)+diff(table_out.SolarAzimuthAngle(i-2:i-1));
        table_out.SZA(i)=table_out.SZA(i-1)+diff(table_out.SZA(i-2:i-1));
        
    else
        
        % bad line is not the last spectrum; average using times before and
        % after
        table_out.Timehhmmss(i)=table_out.Timehhmmss(i-1)+...
            diff([table_out.Timehhmmss(i-1),table_out.Timehhmmss(i+1)])/2;
        
        % SZA, SAA also needs averaging
        table_out.SolarAzimuthAngle(i)=table_out.SolarAzimuthAngle(i-1)+...
            diff([table_out.SolarAzimuthAngle(i-1),table_out.SolarAzimuthAngle(i+1)])/2;
        table_out.SZA(i)=table_out.SZA(i-1)+diff([table_out.SZA(i-1),table_out.SZA(i+1)])/2;
        
    end

end

%% correct other values 
% use time field

% unrelated: fill other NaT from fractionaldate (matlab sometimes fails to
% read valid time fields)
% NaTs for the bad lines have been filled already, any remining NaTs are
% due to matlab messing up, but Fractionalday is fine for those lines
if sum(isnat(table_out.Timehhmmss))>0
    table_out.Timehhmmss(isnat(table_out.Timehhmmss))=...
        ft_to_date(table_out.Fractionalday(isnat(table_out.Timehhmmss))-1,...
        table_out.Year(isnat(table_out.Timehhmmss)));
end

% redo date field
table_out.DateDDMMYYYY=table_out.Timehhmmss;
table_out.DateDDMMYYYY.Format='dd/MM/yyyy';

% redo fd, ft
table_out.Fractionalday(bad_ind)=fracdate(table_out.Timehhmmss(bad_ind))+1;
table_out.Fractionaltime(bad_ind)=...
    (table_out.Fractionalday(bad_ind)-floor(table_out.Fractionalday(bad_ind)))*24;

%%% leave TotalExperimentTimesec as negative, to leave some trace that the
%%% file had bad times


end





