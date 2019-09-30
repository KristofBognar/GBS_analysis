function [ ind_goodvcd, VCD_table ] = filter_VCD_output(tg, VCD_table, rcd_S, instr,...
                                      addRCDerr, do_cams_filter)
%filter_VCD_output(tg,VCD_table,rcd_S,instr,addRCDerr) Filters saved VCD data and returns
%indices of rows with good values
%   INPUT:
%       tg: trace gas selection
%               1 - GBS ozone
%               2,3 - GBS NO2 and NO2 UV
%               4 - SAOZ ozone
%               5 - SAOZ NO2
%       VCD_table, rcd_S -- output of VCD code, saved in a .mat file
%       instr: instrument identifier
%               1 - UT-GBS
%               2 - PEARL-GBS
%               3 - SAOZ with fix RCD
%               4 - SAOZ with daily RCD
%       addRCDerr: if true, code fills in missing GBS systematic errors
%               when the VCD passes all other criteria (one twilight dind't
%               have good RCD fit). If false, days with one twilight are
%               filtered out
%       do_cams_filter: if true, makes sure that chunk of data has less than
%               50% of measurements over/under the CAMS error limits
%   OUTPUT:
%       ind_goodvcd: indices of rows with good values


%% setup

if nargin==4
    addRCDerr=false;
    do_cams_filter=false;
elseif nargin==5
    do_cams_filter=false;
end

% minimum required solar zenith angle range for twilight measurements
sza_range_lim=2;

% min number of points requiredin Langley plot
min_npoints=8;

% extra lines to filter out (selected manually)
% correspond to lines of the unfiltered VCD table
UT_badlines_no2=containers.Map;

% for NO2 small window
% UT_badlines_no2('1999')=[28];
% UT_badlines_no2('2003')=[13:34];
% UT_badlines_no2('2006')=[2,45];
% UT_badlines_no2('2007')=[38];
% UT_badlines_no2('2009')=[3,24,26];
% UT_badlines_no2('2010')=[174];
% UT_badlines_no2('2012')=[134];
% UT_badlines_no2('2013')=[180];
% UT_badlines_no2('2014')=[372];
% UT_badlines_no2('2015')=[195];

% for NO2 large window
UT_badlines_no2('1999')=[28];
UT_badlines_no2('2003')=[21,22];
UT_badlines_no2('2012')=[121:207]; % overheating issues, no operator

% ozone
UT_badlines_o3=containers.Map;
UT_badlines_o3('2005')=[17];
UT_badlines_o3('2007')=[27];
UT_badlines_o3('2009')=[7];
UT_badlines_o3('2010')=[91];
UT_badlines_o3('2011')=[2,453];
UT_badlines_o3('2012')=[131:220]; % overheating issues, no operator
UT_badlines_o3('2013')=[178];
UT_badlines_o3('2014')=[154,404,405];
UT_badlines_o3('2015')=[3];
UT_badlines_o3('2016')=[9];
UT_badlines_o3('2017')=[260,261];

% same for PEARL-GBS
P_badlines_no2=containers.Map;
P_badlines_no2('2008')=[192];

P_badlines_no2uv=containers.Map;
P_badlines_no2uv('2011')=[38,39];
P_badlines_no2uv('2012')=[70:333]; %CCD region was severely restricted, fit quality is bad
P_badlines_no2uv('2013')=[110:205]; % gain set to best dynamic range, data is garbage
P_badlines_no2uv('2017')=[133,134,150,151,211,212,229:232];

P_badlines_o3=containers.Map;
P_badlines_o3('2006')=[70,71];
P_badlines_o3('2008')=[80];
P_badlines_o3('2009')=[25,26];

% same for SAOZ
S_badlines_no2=containers.Map;
S_badlines_no2('2009')=[53:56,72];
S_badlines_no2('2012')=[44];

%% grab info from RCD structure

% check if VCD_table and rcd_S are different length
if length(VCD_table.mean_vcd)~=length(rcd_S.R2)

    % find matching RCD fields -- use fd_min, it should be the same in both
    % VCD_table and rcd_S
    [~,ind_matchvcd,ind_matchrcd]=intersect(VCD_table.fd_min,rcd_S.fd_min);

    % initialize arrays for R2, sza range, and number of points in langley plot
    R2=zeros(size(VCD_table.mean_vcd));
    sza_range=zeros(size(VCD_table.mean_vcd));
    npoints=zeros(size(VCD_table.mean_vcd));

    % part of the year has NaN VCDs: RCDs are only calculated for good
    % VCDs, so assign values accordingly (missing VCDs have R2=0)
    R2(ind_matchvcd)=rcd_S.R2(ind_matchrcd);
    sza_range(ind_matchvcd)=rcd_S.sza_max(ind_matchrcd)-rcd_S.sza_min(ind_matchrcd);
    npoints(ind_matchvcd)=rcd_S.nbr(ind_matchrcd);

else
    R2=rcd_S.R2;
    sza_range=rcd_S.sza_max-rcd_S.sza_min;
    npoints=rcd_S.nbr;
end

% get year info
year=VCD_table.year(1);

%% filter VCDs
if tg==1 || tg==4
    % ozone: filter for NaNs and 60 DU error limit (from Cristen's thesis)

    % set R^2 limit
    % allow more error for 2003-2004, since data is bad quality (bad RMS filter used)
    if year==2003 || year==2004
        R2_lim=0.8;
    elseif instr==3
        R2_lim=0; % exclude R2 limit for fix RCD fit, langley fit not included in retrieval
        sza_range_lim=0.5; % no langley plot required for fix ref data
    else
        R2_lim=0.9;
    end

    % calculate total error
    vcd_error=sqrt(VCD_table.sigma_mean_vcd.^2 + VCD_table.std_vcd.^2);
    
    % good values
    ind_goodvcd=find(~isnan(VCD_table.mean_vcd) & ...
                     VCD_table.mean_vcd>0 & ...
                     vcd_error < 60*2.687e16 & ...
                     R2>=R2_lim & ...
                     sza_range>=sza_range_lim & ...
                     npoints>min_npoints);

    if addRCDerr
        % find good values when the other twilight failed for some reason (sza range,
        % number of points, R2 or simply missing)
        % these VCDs must pass all other error filters
        ind_stillOK=find(~isnan(VCD_table.mean_vcd) & ...
                         VCD_table.mean_vcd>0 & ...
                         isnan(VCD_table.sigma_mean_vcd) & ...
                         ~isnan(VCD_table.std_vcd) & ...
                         R2>=R2_lim & ...
                         sza_range>=sza_range_lim & ...
                         npoints>min_npoints);

    %     ind_tmp=find(isnan(VCD_table.sigma_mean_vcd) & ~isnan(VCD_table.std_vcd));
    %     VCD_table.sigma_mean_vcd(ind_tmp)=7e17;

        % fill missing syst error values: interpolate + extrapolate by keeping
        % values outside the range constant
        VCD_table.sigma_mean_vcd(ind_stillOK)=...
            interp1([1;VCD_table.fd(ind_goodvcd);365],...
                    [VCD_table.sigma_mean_vcd(ind_goodvcd(1));...
                     VCD_table.sigma_mean_vcd(ind_goodvcd);...
                     VCD_table.sigma_mean_vcd(ind_goodvcd(end))],...
                    VCD_table.fd(ind_stillOK));

        % final list of good indices
        ind_goodvcd=sort([ind_stillOK;ind_goodvcd]);
    end
    
    % apply manual filter to outliers
    if instr==1 && isKey(UT_badlines_o3,num2str(year))
        ind_goodvcd=setdiff(ind_goodvcd,UT_badlines_o3(num2str(year)));
    elseif instr==2 && isKey(P_badlines_o3,num2str(year))
        ind_goodvcd=setdiff(ind_goodvcd,P_badlines_o3(num2str(year)));
    end
    
    if do_cams_filter

        % CAMS error limits
        lim_sys=[2.8,10];
% %         lim_rand=[3.5,5];
        lim_rand=[3.5,6];
       
        [ind_goodvcd]=CAMS_filter(ind_goodvcd,VCD_table,lim_sys,lim_rand);
        
    end
    
elseif tg==2 || tg==3 || tg==5
    % no2: filter for NaNs and 2e15 (or 200%) error limit from Cristen's thesis
    % R^2 limit is 0.6 here, and VCDs are only considered if both twilights have
    % data (so NaN errors should be filtered out)
    % same filters for NO2 UV as for NO2 VIS
    % same filter for GBS and SAOZ daily RCD
    if instr~=3
        R2_lim=0.6;
    else % exclude R2 and SZA limit for fix RCD fit -- langley fit not included in retrieval
        R2_lim=0; 
        sza_range_lim=0.5;
    end
        
    % calculate total error
    vcd_error=sqrt(VCD_table.sigma_mean_vcd.^2 + VCD_table.std_vcd.^2);
    
    
    
    ind_goodvcd=find(~isnan(VCD_table.mean_vcd) & ...
                     VCD_table.mean_vcd>0 & ...
                     vcd_error < 2e15 & ...
                     vcd_error./VCD_table.mean_vcd < 2 & ...                     
                     R2>=R2_lim & ...
                     sza_range>=sza_range_lim);

    % apply manual filter to outliers
    if instr==1 && isKey(UT_badlines_no2,num2str(year)) && tg==2
        ind_goodvcd=setdiff(ind_goodvcd,UT_badlines_no2(num2str(year)));
    elseif instr==2 && isKey(P_badlines_no2,num2str(year)) && tg==2
        ind_goodvcd=setdiff(ind_goodvcd,P_badlines_no2(num2str(year)));
    elseif instr==2 && isKey(P_badlines_no2uv,num2str(year)) && tg==3
        ind_goodvcd=setdiff(ind_goodvcd,P_badlines_no2uv(num2str(year)));
    elseif tg==5 && isKey(S_badlines_no2,num2str(year))
        ind_goodvcd=setdiff(ind_goodvcd,S_badlines_no2(num2str(year)));
    end

    if do_cams_filter

        % CAMS error limits
% %         lim_sys=[1.5,24];
% %         lim_rand=[2.8,18];
        lim_sys=[1.5,52];
        lim_rand=[2,40];
       
        [ind_goodvcd]=CAMS_filter(ind_goodvcd,VCD_table,lim_sys,lim_rand);
        
    end
    

    
else
    error('Tracegas input not recognized')
end
 


end

function [ind_goodvcd]=CAMS_filter(ind_goodvcd,VCD_table,lim_sys,lim_rand)

% calculate systematic and random error percentages, and sort them
% in descending order (save indices too)
[sys,ind_sys]=sort((VCD_table.sigma_mean_vcd(ind_goodvcd)...
                   ./VCD_table.mean_vcd(ind_goodvcd))*100,'descend');
[rand,ind_rand]=sort((VCD_table.std_vcd(ind_goodvcd)...
                     ./VCD_table.mean_vcd(ind_goodvcd))*100,'descend');

% see if either sys or rand exceed error limits over 50% of the time
over_lim_sys=sum(sys>=lim_sys(2) | sys<=lim_sys(1));
over_lim_rand=sum(rand>=lim_rand(2) | rand<=lim_rand(1));

lim_50=length(ind_goodvcd)/2;

if over_lim_sys>=lim_50 || over_lim_rand>=lim_50

    % pick error that exceeds limit and n.o. measurements to discard
    if over_lim_sys>over_lim_rand
        ind_over=ind_sys;
        discard=2*over_lim_sys-length(ind_goodvcd)+1;
    else
        ind_over=ind_rand;
        discard=2*over_lim_rand-length(ind_goodvcd)+1;
    end

    if discard>=length(ind_goodvcd)
        ind_goodvcd=[];
        return
    end
    
    % pick out indices to discard
    % ind_sys and ind_rand are indices of ind_goodvcd, and the
    % first <discard> number are the largest errors
    ind_goodvcd(ind_over(1:discard))=[];

else
    return
end
    

end


