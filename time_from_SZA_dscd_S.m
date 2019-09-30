function [ dscd_S_new ] = time_from_SZA_dscd_S( instr, dscd_S )
%time_from_SZA_readQDOAS( instr, dscd_S ) 
%   Finds times of reference spectra from SZA info and adds it to dscd_S structure
%	Timestep is defined in time_from_SZA.m
%
%	Modifications to dscd_S:
%		dscd_S.ref_sza:	-1 values for fixed references replaced using records in manual_refs.mat
%		dscd_S.ref_sza_utc1: One of the possible times for the given SZA (matlab serial time, UTC)
%		dscd_S.ref_sza_utc2: The other possible time for the given SZA (matlab serial time, UTC)
%							 If the first time is <+-30 min from noon (taken as 18:00), then
%							 thi value is NaN and only the first value is reported 						


%% find manual refs
ind_manual=find(dscd_S.ref_sza==-1);

% load manual ref sza file
load manual_refs.mat

% input for time calculation
year=dscd_S.year(ind_manual);
doy=dscd_S.day(ind_manual);

% display details
disp(['Manual references used for ' num2str(length(unique(doy))) ' days'])

% match year and day to manual ref sza
sza_fill=[];
for i=1:length(ind_manual)
    ind=find(manual_refs.instr==instr & manual_refs.year==year(i) & manual_refs.doy==doy(i));
    if ~isempty(ind)
        sza_fill=[sza_fill; manual_refs.sza(ind)];
    else
        % check if all days have consistent ref_SZA values (QDOAS sometimes writes
        % -1 for first spec of next day)
        if dscd_S.ref_sza(ind_manual(i)+1)~=-1 && dscd_S.day(ind_manual(i)+1)==doy(i)
            sza_fill=[sza_fill; dscd_S.ref_sza(ind_manual(i)+1)];
        else
            sza_fill=[sza_fill; 0];
        end
    end
end

% complete dscd_S.ref_sza
dscd_S.ref_sza(ind_manual)=sza_fill;

%% get times from sza values

% find sza for each day
[~,ia,ic] = unique(dscd_S.day);

time=time_from_SZA(dscd_S.year(ia),dscd_S.day(ia),dscd_S.ref_sza(ia),'E');

dscd_S.ref_sza_utc1=time(ic,1);
dscd_S.ref_sza_utc2=time(ic,2);

% return modified structure
dscd_S_new=dscd_S;

end

