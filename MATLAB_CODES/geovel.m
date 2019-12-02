% Computes the horizontal geostropic velocities from interface height 
% with periodic boundary conditions
%
% usage: [ug,vg]=geovel(eta,f,gp,dx,dy)
%
% input  eta: interface height [nx x ny x nz+1]
%        f: Coriolis parameter, g' [nx x ny] 
%        gp: reduced gravity, g' [nz]
%        dx: grid spacing in x [nx x ny] (technically this should be 1/2(x(i+1)-x(i-1)))
%        dy: grid spacing in y [nx x ny]   " 
% output ug: zonal geostrophic velocity [nx x ny x nz]
%        vg: meridional geostrophic velocity [nx x ny x nz]



function [ug,vg]=geovel(eta,f,gp,dx,dy)
  gp3D=permute(repmat(gp',[1,size(eta,1),size(eta,2)]),[2,3,1]);
  f3D=repmat(f,[1,1,length(gp)]);
  M=cumsum(eta(:,:,1:end-1).*gp3D,3);
  [dMdx,dMdy]=grad(M,dx,dy);
  vg=f3D.^(-1).*dMdx;
  ug=-f3D.^(-1).*dMdy;
end



