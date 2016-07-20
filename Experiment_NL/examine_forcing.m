clf
nc = netcdf.open('forcing.nc');

taux = netcdf.getVar(nc, netcdf.inqVarID(nc,'taux') );
taux(taux<-1.e33) = NaN
tauy = netcdf.getVar(nc, netcdf.inqVarID(nc,'tauy') );
tauy(tauy<-1.e33) = NaN

%contourf(taux');
%colorbar()
contourf(tauy');
colorbar()
