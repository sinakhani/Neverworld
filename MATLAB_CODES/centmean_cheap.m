% Square box average that remains centered at all times. When missing values
% occur, the local box size is reduced until all values are well defined.
% Produces NaNs only where NaNs exixted in the original data
% Boundary conditions are treated as doubly-periodic. If you want to avoid periodic BC's,
% just add row and/or column of NaNs at the boundary
%
% usage: Xave=centmean(X,k)
% input  X: data matrix (nxm)
%        k: side length of averaging box - has to be odd
% output Xave: Filtered data matrix (nxm)

function Xave = centmean_cheap(X, k)

kk=(k-1)/2;

% Construct extended data array uing periodic BCs:
Xnew=zeros(size(X,1)+k-1,size(X,2)+k-1);
Xnew(kk+1:end-kk,kk+1:end-kk)=X;
Xnew(1:kk,kk+1:end-kk)=X(end-kk+1:end,:);
Xnew(end-kk+1:end,kk+1:end-kk)=X(1:kk,:);
Xnew(:,1:kk)=Xnew(:,end-2*kk+1:end-kk);
Xnew(:,end-kk+1:end)=Xnew(:,kk+1:2*kk);

% Do the original filtering:
% (All points where convolution includes NaNs will be set to NaN here)
window=ones(2*kk+1,2*kk+1)/(2*kk+1)^2;
Xave=conv2(Xnew,window,'same');    
% for i=kk-1:-1:0                  %%%% Comment this loop for cheap case
%    % successively fill in NaN's using smaller boxes
%    window=ones(2*i+1,2*i+1)/(2*i+1)^2;
%    Xhelp=conv2(Xnew,window,'same');
%    Xave(isnan(Xave))=Xhelp(isnan(Xave));
% end                              %%%% Comment this loop for cheap case

% return array with original shape (remove extensions for periodic BCs)
Xave=Xave(kk+1:end-kk,kk+1:end-kk);
