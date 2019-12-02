% Computes the horizontal gradient of a scalar 2D or 3D field 
% with periodic boundary conditions, using centered differences
%
% usage: [dadx,dady]=grad(a,dx,dy)
%
% input  a: data matrix [nx x ny x nz]
%        dx: grid spacing in x [nx x ny] (technically this should be 1/2(x(i+1)-x(i-1)))
%        dy: grid spacing in y [nx x ny]
% output dadx: x derivative  of a [nx x ny x nz]
%        dady: y derivative  of a [nx x ny x nz]


function [dadx,dady]=grad(a,dx,dy)


% Construct extended data array uing periodic BCs:
anew=zeros(size(a,1)+2,size(a,2)+2,size(a,3));
anew(2:end-1,2:end-1,:)=a;
anew(1,2:end-1,:)=a(end,:,:); anew(end,2:end-1,:)=a(1,:,:);
anew(:,1,:)=anew(:,end-1,:); anew(:,end,:)=anew(:,2,:);

dxnew=zeros(size(a,1)+2,size(a,2)+2);
dxnew(2:end-1,2:end-1)=dx;
dxnew(1,2:end-1)=dx(end,:); dxnew(end,2:end-1)=dx(1,:);
dxnew(:,1)=dxnew(:,end-1); dxnew(:,end)=dxnew(:,2);

dynew=zeros(size(a,1)+2,size(a,2)+2);
dynew(2:end-1,2:end-1)=dy;
dynew(1,2:end-1)=dy(end,:); dynew(end,2:end-1)=dy(1,:);
dynew(:,1)=dynew(:,end-1); dynew(:,end)=dynew(:,2);

dxnew=repmat(dxnew,[1,1,size(a,3)]);
dynew=repmat(dynew,[1,1,size(a,3)]);

dadx=(anew(3:end,2:end-1,:)-anew(1:end-2,2:end-1,:))./(2*dxnew(2:end-1,2:end-1,:));
dady=(anew(2:end-1,3:end,:)-anew(2:end-1,1:end-2,:))./(2*dynew(2:end-1,2:end-1,:));

end



