function [ col ] = create_colnum_struct(str,tg)
% Create structure of column numbers to read data from QDOAS output file
%
% INPUT: string to select appropriate file format ('instr_tg')
%   'utgbs_o3': UT-GBS ozone (same file as NO2)
%   'utgbs_no2': UT-GBS NO2 (same file as ozone)
%   'pgbs_no2': PEARL-GBS NO2 
%   'pgbs_bro': PEARL-GBS MAX-DOAS BrO data
%   'gbs_redo_vis': output from UT-GBS reanalysis QDOAS project (visible range)
%   'gbs_o3x-: O3 data with X.xs from UT-GBS reanalysis project (visible range)
%   'gbs_redo_vis_p': output from PEARL-GBS reanalysis QDOAS project (visible range)
%   'gbs_redo_uv': output from PEARL-GBS AND UT-GBS reanalysis QDOAS project (UV range)
%
%   tg: 1 for ozone, 2 for no2, 3 for o3x; only of gbs_redo_vis is selected
%
% Created for easier modification of the column structure
% Kristof Bognar, June 2017

%% select input option

% create empty structure
col=struct();

if nargin>1 && ~strcmp(str,'gbs_redo_vis') && ~strcmp(str,'gbs_redo_vis_p') && ~strcmp(str,'gbs_redo_uv')
    error('Tracegas selection only relevant for GBS reanalysis files');
end

% initialize selection variables
utgbs_o3=false;
utgbs_no2=false;
pgbs_no2=false;
pgbs_bro=false;
gbs_redo_vis=false;
gbs_redo_vis_p=false;
gbs_redo_uv=false;
gbs_o3x=false;
cristen=false;

redo_version=3;

% regular zenith-sky QDOAS output
if strcmp(str,'utgbs_o3')
    utgbs_o3=true;
elseif strcmp(str,'utgbs_no2')
    utgbs_no2=true;
elseif strcmp(str,'pgbs_no2')
    pgbs_no2=true;
    
% PGBS maxdoas output (BrO only)
elseif strcmp(str,'pgbs_bro')
    pgbs_bro=true;
    
% GBS reanalysis output (visible)
elseif strcmp(str,'gbs_redo_vis')
    gbs_redo_vis=true;

elseif strcmp(str,'gbs_redo_vis_p')
    gbs_redo_vis_p=true;

elseif strcmp(str,'gbs_redo_uv')
    gbs_redo_uv=true;

% O3x output from GBS reanalysis (visible)
elseif strcmp(str,'gbs_o3x')
    gbs_o3x=true;

% Cristen's DSCD files
elseif strcmp(str,'cristen')
    cristen=true;

else error('Incorrect input')
end

%% general fields
if utgbs_no2 || utgbs_o3 || pgbs_no2

    if pgbs_no2, col.tot_nbr=45; else col.tot_nbr=41; end
    
    col.year=2;
    col.fd=3;
    col.ft=4;
    col.sza=7;
    col.saa=8;
    col.elev=9;
    col.tint=5;
    col.scans=6;

    col.ref_sza=13;
    
% elseif gbs_redo || pgbs_bro
elseif pgbs_bro    
    
    if pgbs_bro, col.tot_nbr=65; else col.tot_nbr=263; end
    col.year=2;
    col.fd=3;
    col.ft=4;
    col.sza=5;
    col.saa=6;
    col.elev=7;
    col.tint=0;
    col.scans=0;

elseif gbs_redo_vis || gbs_o3x || cristen || gbs_redo_vis_p || gbs_redo_uv
    
    % reanalysis with QDOAS 2.109 (.L2 files, has tint and scans in output)
    if redo_version==2, col.tot_nbr=198; end

    % reanalysis with QDOAS 3.1 (L2.ASC files, has tint and scans in output)
    if redo_version==3, col.tot_nbr=204; end
    
    if gbs_redo_vis_p, col.tot_nbr=230; end
    
    if gbs_redo_uv, col.tot_nbr=78; end
    
    if gbs_o3x, col.tot_nbr=58; end

    if cristen, col.tot_nbr=42; end

    % same for all
    col.year=2;
    col.fd=3;
    col.ft=4;
    col.sza=7;
    col.saa=8;
    col.elev=9;
    col.tint=6;
    col.scans=5;
    
end

%% tracegas windows

if utgbs_no2

    col.rms=27;
    col.dscd=35;
    col.err=36;
    col.shift=39;
    col.stretch=40;

elseif utgbs_o3
    
    col.rms=12;
    col.dscd=22;
    col.err=23;
    col.shift=24;
    col.stretch=25;
    
elseif pgbs_no2
    
    col.rms=12;
    col.dscd=18;
    col.err=19;
    col.shift=26;
    col.stretch=27;
    col.oclo_dscd=22;
    col.oclo_err=23;
    
elseif pgbs_bro

    col.ref_sza=30;
    
    col.rms=29;
    col.dscd=41;
    col.err=42;
    col.shift=44;
    col.stretch=45;
    col.dscd_o4=49;
    col.err_o4=50;
    
elseif gbs_o3x

    col.ref_sza=15;
    
    col.rms=14;
    col.dscd=25;
    col.err=26;
    col.shift=37;
    col.stretch=38;
    
    col.x=27;

elseif cristen

    col.ref_sza=42;
    
    col.rms=12;
    col.dscd=21;
    col.err=22;
    col.shift=25;
    col.stretch=26;
    
    col.x=23;
    
elseif gbs_redo_vis
    
    if redo_version==3
        % reanalysis with QDOAS 3.1 (L2.ASC files, has tint and scans in output)

        % o4 293K, with 203K orthogonalized
        col.o4_293a203_rms=14;
        col.o4_293a203_dscd=23;
        col.o4_293a203_err=col.o4_293a203_dscd+1;

        % o4 203K, with 293K orthogonalized
        col.o4_203a293_rms=48;
        col.o4_203a293_dscd=57;
        col.o4_203a293_err=col.o4_203a293_dscd+1;

        % o4 203K
        col.o4_203_rms=82;
        col.o4_203_dscd=91;
        col.o4_203_err=col.o4_203_dscd+1;

        % o4 293K
        col.o4_293_rms=114;
        col.o4_293_dscd=123;
        col.o4_293_err=col.o4_293_dscd+1;

        if tg==1
            % o3, with standard 223K c-s
            col.rms=146;
            col.dscd=157;
            col.err=col.dscd+1;

            col.ref_sza=147;
            col.shift=167;
            col.stretch=168;            
            
        elseif tg==2
            % no2, with standard 223K c-s
            col.rms=170;
            col.dscd=179;
            col.err=col.dscd+1;

            col.ref_sza=171;
            col.shift=183;
            col.stretch=184;            
        
        end
        
        % fluxes
        num_f=19;
        cols=num2cell([col.tot_nbr-num_f+1:col.tot_nbr]);
        [col.f355, col.f360, col.f380, col.f385, col.f390, col.f405,...
         col.f420, col.f425, col.f435, col.f440, col.f445,...
         col.f450, col.f455, col.f460, col.f470, col.f490,...
         col.f500, col.f532, col.f550]=deal(cols{:});
 
    elseif redo_version==2
        % reanalysis with QDOAS 2.109 (.L2 files, has tint and scans in output)

        % o4 293K, with 203K orthogonalized
        col.o4_293a203_rms=14;
        col.o4_293a203_dscd=22;
        col.o4_293a203_err=col.o4_293a203_dscd+1;

        % o4 203K, with 293K orthogonalized
        col.o4_203a293_rms=47;
        col.o4_203a293_dscd=55;
        col.o4_203a293_err=col.o4_203a293_dscd+1;

        % o4 203K
        col.o4_203_rms=80;
        col.o4_203_dscd=88;
        col.o4_203_err=col.o4_203_dscd+1;

        % o4 293K
        col.o4_293_rms=111;
        col.o4_293_dscd=119;
        col.o4_293_err=col.o4_293_dscd+1;

        % o3, with standard 223K c-s
        col.o3_rms=142;
        col.o3_dscd=152;
        col.o3_err=col.o3_dscd+1;

        % no2, with standard 223K c-s
        col.no2_rms=165;
        col.no2_dscd=173;
        col.no2_err=col.no2_dscd+1;

        % fluxes
        num_f=19;
        cols=num2cell([col.tot_nbr-num_f+1:col.tot_nbr]);
        [col.f355, col.f360, col.f380, col.f385, col.f390, col.f405,...
         col.f420, col.f425, col.f435, col.f440, col.f445,...
         col.f450, col.f455, col.f460, col.f470, col.f490,...
         col.f500, col.f532, col.f550]=deal(cols{:});
        
    end
    
elseif gbs_redo_vis_p
    
    % reanalysis with QDOAS 3.1 (L2.ASC files, has tint and scans in output)
    % files have both O3 and O3X data

    % o4 293K, with 203K orthogonalized
    col.o4_293a203_rms=14;
    col.o4_293a203_dscd=23;
    col.o4_293a203_err=col.o4_293a203_dscd+1;

    % o4 203K, with 293K orthogonalized
    col.o4_203a293_rms=48;
    col.o4_203a293_dscd=57;
    col.o4_203a293_err=col.o4_203a293_dscd+1;

    % o4 203K
    col.o4_203_rms=82;
    col.o4_203_dscd=91;
    col.o4_203_err=col.o4_203_dscd+1;

    % o4 293K
    col.o4_293_rms=114;
    col.o4_293_dscd=123;
    col.o4_293_err=col.o4_293_dscd+1;

    if tg==1
        % o3, with standard 223K c-s
% o3 with no X.xs        
%         col.rms=146;
%         col.dscd=157;
%         col.err=col.dscd+1;
% 
%         col.ref_sza=147;
%         col.shift=167;
%         col.stretch=168;      
        
        % o3, with standard 223K c-s + X.xs
        col.rms=170;
        col.dscd=181;
        col.err=col.dscd+1;

        col.ref_sza=171;
        col.shift=193;
        col.stretch=194;      

        col.x=183;

    elseif tg==2
        % no2, with standard 223K c-s
        col.rms=196;
        col.dscd=205;
        col.err=col.dscd+1;

        col.ref_sza=197;
        col.shift=209;
        col.stretch=210;            
    end

    % fluxes
    num_f=19;
    cols=num2cell([col.tot_nbr-num_f+1:col.tot_nbr]);
    [col.f355, col.f360, col.f380, col.f385, col.f390, col.f405,...
     col.f420, col.f425, col.f435, col.f440, col.f445,...
     col.f450, col.f455, col.f460, col.f470, col.f490,...
     col.f500, col.f532, col.f550]=deal(cols{:});
    
elseif gbs_redo_uv
    
    % UV reanalysis with QDOAS 3.1 (L2.ASC files, has tint and scans in output)
    % files have NO2 and O3 data

    if tg==1
        % o3, with standard 223K c-s
        col.rms=14;
        col.dscd=21;
        col.err=col.dscd+1;

        col.ref_sza=15;
        col.shift=39;
        col.stretch=40;      

    elseif tg==2
        % no2, with standard 223K c-s
        col.rms=42;
        col.dscd=49;
        col.err=col.dscd+1;

        col.ref_sza=43;
        col.shift=57;
        col.stretch=58;            

    end
    
% elseif gbs_redo
%     
%     
 
% % % elseif gbs_redo
% % %     
% % %     % o4 293K, with 203K orthogonalized
% % %     col.o4_293a203_rms=9;
% % %     col.o4_293a203_dscd=18;
% % %     col.o4_293a203_err=col.o4_293a203_dscd+1;
% % % 
% % %     % o4 203K, with 293K orthogonalized
% % %     col.o4_203a293_rms=43;
% % %     col.o4_203a293_dscd=52;
% % %     col.o4_203a293_err=col.o4_203a293_dscd+1;
% % %     
% % %     % o4 203K
% % %     col.o4_203_rms=77;
% % %     col.o4_203_dscd=86;
% % %     col.o4_203_err=col.o4_203_dscd+1;
% % %     
% % %     % o4 293K
% % %     col.o4_293_rms=109;
% % %     col.o4_293_dscd=118;
% % %     col.o4_293_err=col.o4_293_dscd+1;
% % % 
% % %     % o3, with standard 223K c-s
% % %     col.o3_rms=141;
% % %     col.o3_dscd=152;
% % %     col.o3_err=col.o3_dscd+1;
% % % 
% % %     % o3, with warm 293K c-s
% % %     col.o3_293_rms=165;
% % %     col.o3_293_dscd=176;
% % %     col.o3_293_err=col.o3_293_dscd+1;
% % % 
% % %     % o3, with offline ring calculation (and strandard 223K c-s)
% % %     col.o3_offline_rms=189;
% % %     col.o3_offline_dscd=200;
% % %     col.o3_offline_err=col.o3_offline_dscd+1;
% % % 
% % %     % no2, with standard 223K c-s
% % %     col.no2_rms=213;
% % %     col.no2_dscd=222;
% % %     col.no2_err=col.no2_dscd+1;
% % % 
% % %     % no2, with warm 298K c-s
% % %     col.no2_298_rms=229;
% % %     col.no2_298_dscd=238;
% % %     col.no2_298_err=col.no2_298_dscd+1;
% % % 
% % %     % fluxes
% % %     num_f=19;
% % %     cols=num2cell([col.tot_nbr-num_f+1:col.tot_nbr]);
% % %     [col.f340, col.f347, col.f350, col.f355, col.f360, col.f380, col.f385,...
% % %      col.f390, col.f405, col.f420, col.f425, col.f435, col.f440, col.f445,...
% % %      col.f450, col.f490, col.f500, col.f532, col.f550]=deal(cols{:}); 
end
