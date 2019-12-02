% Function to compute SGS volume transport from interface height 
% with periodic boundary conditions
%
% usage: [ugh_SGS,vgh_SGS]=SGSflux(eta,f,gp,dx,dy,kfilt,eps,crops)
%
% input  eta: interface height [nx x ny x nz+1]
%         f: Coriolis parameter, g' [nx x ny] 
%        gp: reduced gravity, g' [nz]
%        dx: grid spacing in x [nx x ny] (technically this should be 1/2(x(i+1)-x(i-1)))
%        dy: grid spacing in y [nx x ny]   " 
%     kfilt: filter length - in gridpoints
%       eps: minimum thickness below which layer is assumed to be vanished [default=0]
%     crops: Treatment of isopycnal in- and out-crops. Options:
%            0 [default] include both surface and bottom in/out-crops as physical layers with zero thickness
%            1 include only surface outcrops as physical layers with zero thickness
%            notice that treating outcrops as boundaries doesn't make 
%            as much sense here since all geostrophic velocities below
%            would become ill defined. (Unless we are somehow treating free surface separately) 
%            
% output ug: zonal geostrophic velocity [nx x ny x nz]
%        vg: meridional geostrophic velocity [nx x ny x nz]

function [ugh_SGS,vgh_SGS,h_ave,ugh_fxy,ug_fxyh_fxy,vgh_fxy,vg_fxyh_fxy]=SGSflux_filtfirst_cheap(eta,f,gp,dx,dy,kfilt,eps,crops)

if nargin<8
    crops=0;
    eps=0; % notice that eps is not actually used in this case
end

eta(eta<-1e33)=NaN;

h=eta(:,:,1:end-1)-eta(:,:,2:end);

Mask=ones(size(h)); 
Mask(isnan(eta(:,:,1:end-1)))=NaN;
Mask_eta=ones(size(eta));Mask_eta(:,:,1:end-1)=Mask;
if crops==1
  % mask out in-cropped isopycnals:
  Mask(eta(:,:,1:end-1)-repmat(eta(:,:,end),[1,1,size(h,3)])<eps)=NaN;
  Mask_eta(:,:,1:end-1)=Mask;
elseif crops ~= 0
   error('crops has to be 0 or 1') 
end

[ug,vg]=geovel(Mask_eta.*eta,f,gp,dx,dy);

A=dx.*dy;

eta_ave=eta; ugh_ave=ug;vgh_ave=vg;
for k=1:size(h,3)
 A_ave=centmean_cheap(squeeze(Mask(:,:,k).*A(:,:)),kfilt);
 ugh_ave(:,:,k)=centmean_cheap(squeeze(Mask(:,:,k).*ug(:,:,k).*h(:,:,k)).*A,kfilt)./A_ave;
 vgh_ave(:,:,k)=centmean_cheap(squeeze(Mask(:,:,k).*vg(:,:,k).*h(:,:,k)).*A,kfilt)./A_ave;
 eta_ave(:,:,k)=centmean_cheap(squeeze(Mask_eta(:,:,k).*eta(:,:,k)).*A,kfilt)./A_ave;
end
A_ave=centmean_cheap(squeeze(Mask_eta(:,:,end).*A(:,:)),kfilt);
eta_ave(:,:,end)=centmean_cheap(squeeze(Mask_eta(:,:,end).*eta(:,:,end)).*A,kfilt)./A_ave;

h_ave=eta_ave(:,:,1:end-1)-eta_ave(:,:,2:end);

[ug_ave,vg_ave]=geovel(eta_ave,f,gp,dx,dy);
ug_ave=Mask.*ug_ave; vg_ave=Mask.*vg_ave;

ugh_SGS=ugh_ave-ug_ave.*h_ave;
vgh_SGS=vgh_ave-vg_ave.*h_ave;
ugh_fxy = ugh_ave;
vgh_fxy = vgh_ave;
ug_fxyh_fxy = ug_ave.*h_ave;
vg_fxyh_fxy = vg_ave.*h_ave;

end
