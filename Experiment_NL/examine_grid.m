nc = netcdf.open('ocean_geometry.nc');

subplot(221)
varid = netcdf.inqVarID(nc,'D');
D = netcdf.getVar(nc,varid);
contourf(D'); colorbar();

subplot(222)
dx = netcdf.getVar(nc, netcdf.inqVarID(nc,'dxT') );
dy = netcdf.getVar(nc, netcdf.inqVarID(nc,'dyT') );
contourf(dx'); colorbar();

subplot(223)
contourf(dy'); colorbar();

subplot(224)
contourf(dx'./dy'-1.); colorbar();

latq = netcdf.getVar(nc, netcdf.inqVarID(nc,'latq') );
lath = netcdf.getVar(nc, netcdf.inqVarID(nc,'lath') );
display( [min(latq), max(latq), min(lath), max(lath)])
lonq = netcdf.getVar(nc, netcdf.inqVarID(nc,'lonq') );
lonh = netcdf.getVar(nc, netcdf.inqVarID(nc,'lonh') );
display( [min(lonq), max(lonq), min(lonh), max(lonh)])
