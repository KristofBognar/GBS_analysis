function BL_out = get_BL_height( year )
% Calculate boundary layer height from sonce PTU data

BL_out=table();

if size(year,1)~=1, year=year'; end

for i=1:length(year)
    
    tmp=table();
    
    % load P, T data for given year
    load(['/home/kristof/work/radiosonde/Eureka/radiosonde_' num2str(year(i)) '_interp.mat']);

    % launch times indicated as noon and midnight in GRAW files; 
    % launches are probably 11 someting and 23 something -- subtract 30 min
    % from fractional times
    if year(i)>=2019; x=x-0.01; end

    pt_times=ft_to_date(x,year(i)); 

    pt_alt=y;
    
    % load wind data for given year
    load(['/home/kristof/work/radiosonde/Eureka/radiosonde_wnd_' num2str(year(i)) '_interp.mat']);

    if year(i)>=2019; x=x-0.01; end

    wnd_times=ft_to_date(x,year(i)); 
    
    wnd_alt=y;
    
    
    % get common times, in case data files are corrupted
    [times,pt_ind,wnd_ind]=intersect(pt_times,wnd_times);

    % remove profiles with only wind or only P, T data
    P_arr=P_arr(:,pt_ind);
    T_arr=T_arr(:,pt_ind);
    wdir_arr=wdir_arr(:,wnd_ind);
    wspd_arr=wspd_arr(:,wnd_ind);
    
    % calculate derived values
    
    % perpendicular wind directions
    wnd_u=wspd_arr.*sind(wdir_arr);
    wnd_v=wspd_arr.*cosd(wdir_arr);
    
    % potential temperature
    theta=get_theta(P_arr, T_arr);

    % get BL height
    bl_height=richardson_BL(pt_alt, theta, wnd_u, wnd_v);
    
    tmp.DateTime=times';
    tmp.BL_height_m=bl_height';
    
    BL_out=[BL_out;tmp];
    
end

end

function bl_height=richardson_BL(alt, theta, wnd_u, wnd_v, Ri_lim, z_s)
% estimate BL height using Eqn 2 of Vogelezang and Holtslag, 1996
% this is thre Bulk Richardson Number for a layer that starts near
% the surface and ends at the top of the BL
    
    bl_height=NaN(1,size(theta,2));

    % defaults
    g=9.81;
    if nargin==4

        Ri_lim=0.3; % works best across tested conditions in Vogelezang and Holtslag, 1996
        z_s=15; % not too important, paper tested 20,40,80 meters with
        % similar results
        ind_z_s=1;
        
    else
        
        ind_z_s=find(alt==z_s);
        if isempty(ind_z_s), error('implement interpolation'), end
        
    end
    
    % 'surface' layer values, repeated for each altitude layer above z_s
    theta_s=repmat(theta(ind_z_s,:),length(alt)-ind_z_s,1);
    wnd_u_s=repmat(wnd_u(ind_z_s,:),length(alt)-ind_z_s,1);
    wnd_v_s=repmat(wnd_v(ind_z_s,:),length(alt)-ind_z_s,1);
    
    % all potential BL top layer values (h>z_s)
    theta_h=theta(ind_z_s+1:end,:);
    wnd_u_h=wnd_u(ind_z_s+1:end,:);
    wnd_v_h=wnd_v(ind_z_s+1:end,:);
    
    
    % Richardson number
    Ri=( (g./theta_s).*(theta_h-theta_s).*repmat(alt(ind_z_s+1:end)-z_s,1,size(theta,2)) ) ./...
       ( (wnd_u_h - wnd_u_s).^2 + (wnd_v_h - wnd_v_s).^2 );
    
    Ri(isinf(Ri))=NaN;
    Ri(Ri>Ri_lim)=NaN;
    
    for i=1:size(theta,2)
        
        tmp=find(isnan(Ri(:,i)));
        
        bl_height(i)=alt(ind_z_s+tmp(1)-1);
        
    end
   
end






