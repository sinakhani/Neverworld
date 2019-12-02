% Function to compute masked geostrophic velocity for spectrum 
% with periodic boundary conditions
%
% input  eta: interface height [nx x ny x nz+1]
%         f: Coriolis parameter, g' [nx x ny] 
%        gp: reduced gravity, g' [nz]
%        dx: grid spacing in x [nx x ny] (technically this should be 1/2(x(i+1)-x(i-1)))
%        dy: grid spacing in y [nx x ny]   " 
%       eps: minimum thickness below which layer is assumed to be vanished [default=0]
%     crops: Treatment of isopycnal in- and out-crops. Options:
%            0 [default] include both surface and bottom in/out-crops as physical layers with zero thickness
%            1 include only surface outcrops as physical layers with zero thickness
%            2 treat all in- or out-crops as boundaries
%            
% output ug: zonal geostrophic velocity [nx x ny x nz]
%        vg: meridional geostrophic velocity [nx x ny x nz]

function [ug_m,vg_m]=spc_geovel(eta,f,gp,dx,dy,eps,crops)

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

ug_m = ug;
vg_m = vg;
end
