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
%            2 treat all in- or out-crops as boundaries
%            
% output ug: zonal geostrophic velocity [nx x ny x nz]
%        vg: meridional geostrophic velocity [nx x ny x nz]

function [ugug_SGS,vgvg_SGS]=SGSeddvis_cheap(eta,f,gp,dx,dy,kfilt,eps,crops)

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
elseif crops==2
  % Mask out all isopycnal in- and out-crops in filter mask:
  Mask(h<eps)=NaN;
  % Mask_eta is not used for filter but only to avoid finite differences
  % into topography in ug,vg calculation. Hence we only need to
  % eliminate incrops here. (Setting surface outcrops to NaN here
  % would lead to all geostrophic velocities below being NaN.)
  helpMask=ones(size(h));
  helpMask(eta(:,:,1:end-1)-repmat(eta(:,:,end),[1,1,size(h,3)])<eps)=NaN;
  Mask_eta(:,:,1:end-1)=helpMask;
elseif crops ~= 0
    error('crops needs to be 0, 1, or 2')
end

[ug,vg]=geovel(Mask_eta.*eta,f,gp,dx,dy);

A=dx.*dy;

ug_ave=ug;vg_ave=vg; ugug_ave=ug;vgvg_ave=vg;
for k=1:size(h,3)
 A_ave=centmean_cheap(squeeze(Mask(:,:,k).*A(:,:)),kfilt);
 ug_ave(:,:,k)=centmean_cheap(squeeze(Mask(:,:,k).*ug(:,:,k)).*A,kfilt)./A_ave;
 vg_ave(:,:,k)=centmean_cheap(squeeze(Mask(:,:,k).*vg(:,:,k)).*A,kfilt)./A_ave;
 ugug_ave(:,:,k)=centmean_cheap(squeeze(Mask(:,:,k).*ug(:,:,k).*ug(:,:,k)).*A,kfilt)./A_ave;
 vgvg_ave(:,:,k)=centmean_cheap(squeeze(Mask(:,:,k).*vg(:,:,k).*vg(:,:,k)).*A,kfilt)./A_ave;
end

ugug_SGS=ugug_ave-ug_ave.*ug_ave;
vgvg_SGS=vgvg_ave-vg_ave.*vg_ave;

end
