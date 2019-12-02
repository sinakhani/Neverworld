%%%%%%%% THIS SCRIPT GENERATES SUBGRID SCALE EDDY FLUXES IN NEVERWORLD 
%%%%%%%% CONFIGRATIONS USING MOM6 NETCDF DATA. PLEASE SEE THE PAPER:
%%%%%%%% (https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019MS001721)
%%%%%%%% THIS SCRIPT IS WRITTEN BY SINA KHANI (skhani@uw.edu).
%%%%%%%% THIS CODE IS FREE TO USE SUBJECT TO PROPER CITATION OF THE WORK.

clear all; 
openstr = 'ocean_geometry.nc';  %%%%% READ OCEAN GEOMETRY NETCDF FILE HERE.
fi = netcdf.open(openstr, 'nowrite');
% get the varibles to plot, there are mavny to choose from
% this for loop will tell you where each variable is stored,
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(fi);
for k=0:(numvars-1)
[vn xt dimid nat] = netcdf.inqVar(fi,k);
disp(sprintf('%d: %s',k,vn));
end

DX=netcdf.getVar(fi,10);
DY = netcdf.getVar(fi,11);
AR = DX.* DY;
ARM = mean(AR,1)';


openstr = 'snapshots.nc'; %%%% READ SNAPSHOT NETCDF FILE HERE.
fi = netcdf.open(openstr, 'nowrite');
% get the varibles to plot, there are mavny to choose from
% this for loop will tell you where each variable is stored,
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(fi);
for k=0:(numvars-1)
[vn xt dimid nat] = netcdf.inqVar(fi,k);
disp(sprintf('%d: %s',k,vn));
end

SZ = size(DX);
xsize = SZ(1);
ysize = SZ(2);
zsize = 6;

% RD1=netcdf.getVar(fi,20);
% RD1(RD1<-1.e33) = NaN;
% RD1_ta = nanmean(RD1,3);
% RD1_taza1 = nanmean(RD1_ta,1);
% RD1_taza = permute(RD1_taza1,[2,1]);



tsize = 1825; %%730; 1825%% Maximum number of years for <= 1/4 deg tsize=3650 and for 1/16 tsize=730, tsize=50 for 1/8 deg. 730 for 1/16 %%% 1825 for 1/8
nyr = 10;   %% Do not change this at all! 
avyr = 0.25;  %% !!!!!! for 1/8 in movmean use 0.25 and for 1/16 and 1/4 (3650) use 0.25/4 
TS = 1460; %%1460;%%365;%floor( tsize * ((nyr-avyr)/nyr) ); %% time average for TS<t<tsize

rho_0 = 1000;
grav  = 9.81;
gprime = [9.81,0.015,0.005,0.003,0.002,0.001];
%gprime = gprime'; 
denst = [1028.1,1028.0,1027.8,1027.5,1027.0,1025.5];
rho_st =[1025.5,1027.0,1027.5,1027.8,1028.0,1028.1];
OMEGA = 7.2921e-5;

AR_6 = zeros(xsize,ysize,zsize,'single');

for p=1:zsize
    AR_6 (:,:,p) = AR(:,:);
end

tic

X=netcdf.getVar(fi,0);
Y=netcdf.getVar(fi,1);
Z=netcdf.getVar(fi,2);
CORCOEF = zeros(xsize,ysize,1);
CORCOEF1 = zeros(ysize,1);

y = Y * (111.32 * 1000); %[m]
x = X * (111.32 * 1000); % [m]

z = [0 -150 -400 -1000 -2000 -3000];
zz = fliplr(z);
DeltaZ = [-150,-250,-600,-1000,-1000,-1000];
DeltaZ = DeltaZ';

N2= - gprime' ./ DeltaZ;
Height =[150, 400, 1000, 2000, 3000,4000]; 
Height = Height';
for j=1:ysize
    CORCOEF(:,j) = 2 * OMEGA * sind(Y(j));  %2 * OMEGA * sind(Y(j)); % -4.985E-5 for Y=-20 % -1.1172E-4 for Y=-50;
    CORCOEF1(j) = 2 * OMEGA * sind(Y(j));
end
DEFORL_min = min(sqrt(N2).*Height / max(CORCOEF1));
DEFORL_max = max(sqrt(N2).*Height / min(CORCOEF1));

SLOPE_isn_max = max(abs(DeltaZ)) / min(DY(:));
SLOPE_isn_min = min(abs(DeltaZ)) / max(DY(:));


nfx = 31;
nfy = 31;
nfxy= nfx;
xsize_f = xsize;
ysize_f = ysize;

xsize_pr = xsize+nfx;

beta = zeros(ysize_f,1,'single');

u = zeros(xsize,ysize,zsize,1,'single');
v = zeros(xsize,ysize,zsize,1,'single');
h = zeros(xsize,ysize,zsize,1,'single');
e = zeros(xsize,ysize,zsize,1,'single');

x_f     = zeros(xsize_f,1);
y_f     = zeros(ysize_f,1);
CORCOEF_f = zeros(ysize_f,1);



e_fxy = zeros(xsize_f,ysize_f,zsize,'single');
vg_fxy = zeros(xsize_f,ysize_f,zsize,'single');
ug_fxy = zeros(xsize_f,ysize_f,zsize,'single');
h_fxy = zeros(xsize_f,ysize_f,zsize,'single');
u_fxy = zeros(xsize_f,ysize_f,zsize,'single');
v_fxy = zeros(xsize_f,ysize_f,zsize,'single');
vgh_fxy = zeros(xsize_f,ysize_f,zsize,'single');
ugh_fxy = zeros(xsize_f,ysize_f,zsize,'single');
ugh_fxyh_fxy = zeros(xsize_f,ysize_f,zsize,'single');
KHdyh_fxy = zeros(xsize_f,ysize_f,zsize,'single');
KHdxh_fxy = zeros(xsize_f,ysize_f,zsize,'single');
AR6_fxy  = zeros(xsize_f,ysize_f,zsize,'single');



vgh_sgs = zeros(xsize_f,ysize_f,zsize,'single');
ugh_sgs = zeros(xsize_f,ysize_f,zsize,'single');
vgh_sgsGM = zeros(xsize_f,ysize_f,zsize,'single');
ugh_sgsGM = zeros(xsize_f,ysize_f,zsize,'single');
vghb_sgs = zeros(xsize_f,ysize_f,zsize,'single');
ughb_sgs = zeros(xsize_f,ysize_f,zsize,'single');
dedx = zeros(xsize_f,ysize_f,zsize,'single');
dhdx = zeros(xsize_f,ysize_f,zsize,'single');
dedy = zeros(xsize_f,ysize_f,zsize,'single');
dhdy = zeros(xsize_f,ysize_f,zsize,'single');
dhdyMO = zeros(xsize_f,ysize_f,zsize,'single');
lhdhdy = zeros(xsize_f,ysize_f,zsize,'single');
lhdhdyMO = zeros(xsize_f,ysize_f,zsize,'single');
dvdx = zeros(xsize_f,ysize_f,zsize,'single');
dvdy = zeros(xsize_f,ysize_f,zsize,'single');
DEDX2 = zeros(xsize_f,ysize_f,zsize,'single');
DEDY2 = zeros(xsize_f,ysize_f,zsize,'single');
DHDX2 = zeros(xsize_f,ysize_f,zsize,'single');
DHDY2 = zeros(xsize_f,ysize_f,zsize,'single');
d2edx2 = zeros(xsize_f,ysize_f,zsize,'single');
d2hdx2 = zeros(xsize_f,ysize_f,zsize,'single');
d2edy2 = zeros(xsize_f,ysize_f,zsize,'single');
d2hdy2 = zeros(xsize_f,ysize_f,zsize,'single');
d2hdy2MO = zeros(xsize_f,ysize_f,zsize,'single');
DX_fxy = zeros(xsize_f,ysize_f,'single');
DY_fxy = zeros(xsize_f,ysize_f,'single');
%mask_fxy = zeros(xsize_f,ysize_f,zsize,'single');

Pressta = zeros(xsize,ysize,zsize,'single');
vgta = zeros(xsize,ysize,zsize,'single');
ugta = zeros(xsize,ysize,zsize,'single');
hta = zeros(xsize,ysize,zsize,'single');
eeta = zeros(xsize,ysize,zsize+1,'single');
eta_fxy = zeros(xsize_f,ysize_f,zsize,'single');
ugta_fxy = zeros(xsize_f,ysize_f,zsize,'single');
vgta_fxy = zeros(xsize_f,ysize_f,zsize,'single');
uta_fxy = zeros(xsize_f,ysize_f,zsize,'single');
vta_fxy = zeros(xsize_f,ysize_f,zsize,'single');
hta_fxy = zeros(xsize_f,ysize_f,zsize,'single');
vgh_fxyta = zeros(xsize_f,ysize_f,zsize,'single');
ugh_fxyta = zeros(xsize_f,ysize_f,zsize,'single');
vg_fxyh_fxyta= zeros(xsize_f,ysize_f,zsize,'single');
ug_fxyh_fxyta = zeros(xsize_f,ysize_f,zsize,'single');
vgh_sgsta  = zeros(xsize_f,ysize_f,zsize,'single');
ugh_sgsta  = zeros(xsize_f,ysize_f,zsize,'single');
vgh_sgs0ta  = zeros(xsize_f,ysize_f,zsize,'single');
ugh_sgs0ta  = zeros(xsize_f,ysize_f,zsize,'single');
vgh_sgs2ta  = zeros(xsize_f,ysize_f,zsize,'single');
ugh_sgs2ta  = zeros(xsize_f,ysize_f,zsize,'single');
vgh_sgsGMta = zeros(xsize_f,ysize_f,zsize,'single');
ugh_sgsGMta = zeros(xsize_f,ysize_f,zsize,'single');
ugh_sgs_fluxta = zeros(xsize_f,ysize_f,zsize,'single');
vgh_sgs_fluxta = zeros(xsize_f,ysize_f,zsize,'single');
vghb_sgsta  = zeros(xsize_f,ysize_f,zsize,'single');
ughb_sgsta  = zeros(xsize_f,ysize_f,zsize,'single');
vghb_sgs_fluxta = zeros(xsize_f,ysize_f,zsize,'single');
vghbp_sgs_fluxta  = zeros(xsize_f,ysize_f,zsize,'single');
ughb_sgs_fluxta = zeros(xsize_f,ysize_f,zsize,'single');
nom_sgsta = zeros(xsize_f,ysize_f,zsize,'single');
denom_sgsta = zeros(xsize_f,ysize_f,zsize,'single');
kappa_sgsta = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsta = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsT1ta = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsT2ta = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOT2ta = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOta = zeros(xsize_f,ysize_f,zsize,'single');
nombS_sgsta = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsta = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsT1ta = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsT2ta = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsMOT2ta = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsMOta = zeros(xsize_f,ysize_f,zsize,'single');
Kappa_edta   = zeros(xsize_f,ysize_f,zsize,'single');
Kappa_edMOta = zeros(xsize_f,ysize_f,zsize,'single');
ddx_ughsgsta = zeros(xsize_f,ysize_f,zsize,'single');
ddy_vghsgsta = zeros(xsize_f,ysize_f,zsize,'single');
ddy_ughsgsta = zeros(xsize_f,ysize_f,zsize,'single');
ddx_vghsgsta = zeros(xsize_f,ysize_f,zsize,'single');
ddx_ughsgs_fluxta = zeros(xsize_f,ysize_f,zsize,'single');
ddy_vghsgs_fluxta = zeros(xsize_f,ysize_f,zsize,'single');
dhdyta = zeros(xsize_f,ysize_f,zsize,'single');
dhdyMOta = zeros(xsize_f,ysize_f,zsize,'single');
lhdhdyta = zeros(xsize_f,ysize_f,zsize,'single');
lhdhdyMOta = zeros(xsize_f,ysize_f,zsize,'single');
dedyta = zeros(xsize_f,ysize_f,zsize,'single');
dhdxta = zeros(xsize_f,ysize_f,zsize,'single');
dedxta = zeros(xsize_f,ysize_f,zsize,'single');
d2hdy2ta = zeros(xsize_f,ysize_f,zsize,'single');
d2hdy2MOta = zeros(xsize_f,ysize_f,zsize,'single');
d2edy2ta = zeros(xsize_f,ysize_f,zsize,'single');
d2hdx2ta = zeros(xsize_f,ysize_f,zsize,'single');
d2edx2ta = zeros(xsize_f,ysize_f,zsize,'single');
vgtahta_fxy = zeros(xsize_f,ysize_f,zsize,'single');
ugtahta_fxy = zeros(xsize_f,ysize_f,zsize,'single');
DHDY2ta = zeros(xsize_f,ysize_f,zsize,'single');
DHDY2MOta = zeros(xsize_f,ysize_f,zsize,'single');
DHDX2ta = zeros(xsize_f,ysize_f,zsize,'single');
vgvg_sgsta  = zeros(xsize_f,ysize_f,zsize,'single');
ugug_sgsta  = zeros(xsize_f,ysize_f,zsize,'single');
sum_divta = zeros(xsize_f,ysize_f,zsize,'single');
sum_diff_hta = zeros(xsize_f,ysize_f,zsize,'single');
sum_diff_hMOta = zeros(xsize_f,ysize_f,zsize,'single');
sum_diff_eta = zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsT2ta = zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsMOT2ta = zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsT2ta = zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsMOT2ta = zeros(xsize_f,ysize_f,zsize,'single');
VEL_sgsta = zeros(xsize_f,ysize_f,zsize,'single');

nom_sgsta_G = zeros(xsize_f,ysize_f,zsize,'single');
denom_sgsta_G = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsta_G = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsT1ta_G = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsT2ta_G = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsT2ta_TRG = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsT2ta_STG = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOT2ta_G = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOT2ta_TRG = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOT2ta_STG = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsta_TRG = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsta_STG = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOta_G = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOta_TRG = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOta_STG = zeros(xsize_f,ysize_f,zsize,'single');
nombS_sgsta_G = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsta_G = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsT1ta_G = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsT2ta_G = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsMOT2ta_G = zeros(xsize_f,ysize_f,zsize,'single');
denomh_sgsMOta_G = zeros(xsize_f,ysize_f,zsize,'single');
Kappa_edta_G = zeros(xsize_f,ysize_f,zsize,'single');
Kappa_edMOta_G = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsT2_TRtaG = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOT2_TRtaG = zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsT2ta_G = zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsT2_TRtaG= zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsT2ta_STG  = zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsMOT2ta_G = zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsMOT2_TRtaG= zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsMOT2ta_STG  = zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsT2ta_G = zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsT2_TRtaG= zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsT2ta_STG  = zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsMOT2ta_G = zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsMOT2_TRtaG= zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsMOT2ta_STG  = zeros(xsize_f,ysize_f,zsize,'single');
VEL_sgsta_TRG  = zeros(xsize_f,ysize_f,zsize,'single');
VEL_sgsta_STG  = zeros(xsize_f,ysize_f,zsize,'single');
VEL_sgsta_G    = zeros(xsize_f,ysize_f,zsize,'single');

sea_surf_height = ncread(openstr,'e',[1,1,1,3],[xsize,ysize,1,396]); %% 396
bathymetry      = ncread(openstr,'e',[1,1,7,3],[xsize,ysize,1,1]); %% 396 

sea_surf_height(sea_surf_height<-1.e33) = NaN;
sea_surf_height(isnan(sea_surf_height))=0;
bathymetry(bathymetry<-1.e33) = NaN;
bathymetry(isnan(bathymetry))=0;

MASK_vghsgs = ones(xsize_f,ysize_f,zsize,'single');
MASK_hfxy = ones(xsize_f,ysize_f,zsize,'single');

tic
eps=5;
sum_av = 0;

for i= 1:tsize
        if (i>=TS)
            u = ncread(openstr,'u',[1,1,1,i],[xsize,ysize,zsize,1]);
            v = ncread(openstr,'v',[1,1,1,i],[xsize,ysize,zsize,1]);
            h = ncread(openstr,'h',[1,1,1,i],[xsize,ysize,zsize,1]);
            e = ncread(openstr,'e',[1,1,1,i],[xsize,ysize,zsize,1]);
            ee = ncread(openstr,'e',[1,1,1,i],[xsize,ysize,zsize+1,1]);
            eee = ncread(openstr,'e',[1,1,1,i],[xsize,ysize,zsize+1,1]);
            
            
            u(u<-1.e33) = NaN;
            v(v<-1.e33) = NaN;
            h(h<-1.e33) = NaN;
            e(e<-1.e33) = NaN;
            ee(ee<-1.e33) = NaN;
            eee(eee<-1.e33) = NaN;
            
            
            %%%%%% horiznotal direction (x-y) filtering ---- START!!!!  
            x_f   = x; 
            y_f   = movmean(y,nfy,1);
            CORCOEF_f = movmean(CORCOEF1,nfy,1);
                                         
            for m=1:zsize
                AR6_fxy = centmean_cheap(AR_6(:,:,m),nfxy);
                e_fxy(:,:,m) = centmean(e(:,:,m) .* AR_6(:,:,m),nfxy) ./ AR6_fxy;
            end
           
            DX_fxy = movmean(DX,nfy,2);
            DY_fxy = movmean(DY,nfy,2); 
            AR_fxy = movmean(AR,nfy,2); 
            ARM_fxy = movmean(ARM,nfy,1); 
            %%%%%% horiznotal direction (x-y) filtering ---- END!!!!
            DYY_fxy = mean(DY_fxy,1)';
            DXX_fxy = mean(DX_fxy,1)';
 
            %%%%%%%%%%%%%%%%%%% SET NANs to zero          
            h(isnan(h))=0;
            u(isnan(u))=0;
            v(isnan(v))=0;
            e(isnan(e))=0;
            %%%%%%%%%%%%%%%%%%%
            
            for m=1:ysize_f
                if (m==1)
                    beta(m) = (CORCOEF_f(m+1) - CORCOEF_f(m)) / (DYY_fxy(m));
                elseif (m>1) && (m<ysize)
                    beta(m) = (CORCOEF_f(m+1) - CORCOEF_f(m-1)) / (2*DYY_fxy(m));
                elseif (m==ysize)
                    beta(m) = (CORCOEF_f(m) - CORCOEF_f(m-1)) / (DYY_fxy(m)); 
                end    
            end
            
            [ugh_sgs,vgh_sgs,h_fxy,ugh_fxy,ug_fxyh_fxy,vgh_fxy,vg_fxyh_fxy]=SGSflux_cheap_NV(eee,CORCOEF,gprime,DX,DY,nfxy,eps,1); %%% crops 1 is the best option 
            [ugh_sgs0,vgh_sgs0,h_fxy0,ugh_fxy0,ug_fxyh_fxy0,vgh_fxy0,vg_fxyh_fxy0]=SGSflux_cheap_NV(eee,CORCOEF,gprime,DX,DY,nfxy,eps,0);
            [ugh_sgs2,vgh_sgs2,h_fxy2,ugh_fxy2,ug_fxyh_fxy2,vgh_fxy2,vg_fxyh_fxy2]=SGSflux_cheap_NV(eee,CORCOEF,gprime,DX,DY,nfxy,eps,0);
            
            [ugug_sgs,vgvg_sgs]=SGSeddvis_cheap_NV(eee,CORCOEF,gprime,DX,DY,nfxy,eps,1);
            
            MASK_vghsgs(isnan(vgh_sgs))=0;
            MASK_hfxy(isnan(h_fxy))=0;
            %MASK_reg  = MASK_vghsgs .* MASK_hfxy;        
            
            ugh_sgs(isnan(ugh_sgs))=0;
            vgh_sgs(isnan(vgh_sgs))=0;             
            ugh_sgs0(isnan(ugh_sgs0))=0;
            vgh_sgs0(isnan(vgh_sgs0))=0;             
            ugh_sgs2(isnan(ugh_sgs2))=0;
            vgh_sgs2(isnan(vgh_sgs2))=0;
            h_fxy(isnan(h_fxy))=0;
            ugh_fxy(isnan(ugh_fxy))=0;
            ug_fxyh_fxy(isnan(ug_fxyh_fxy))=0;
            vgh_fxy(isnan(vgh_fxy))=0;
            vg_fxyh_fxy(isnan(vg_fxyh_fxy))=0;
                        

            ugug_sgs(isnan(ugug_sgs))=0;
            vgvg_sgs(isnan(vgvg_sgs))=0;
            
            EKE_sgs = ( ugug_sgs + vgvg_sgs );
            VEL_sgs = real(EKE_sgs.^0.5);
            

                        
            vgh_sgs_flux(:,:,1) =  vgh_sgs(:,:,6);
            ugh_sgs_flux(:,:,1) =  ugh_sgs(:,:,6);

            for k=2:zsize
                vgh_sgs_flux(:,:,k) =  ( vgh_sgs(:,:,zsize+1-k) ) + vgh_sgs_flux(:,:,k-1);
                ugh_sgs_flux(:,:,k) =  ( ugh_sgs(:,:,zsize+1-k) ) + ugh_sgs_flux(:,:,k-1);    
            end
            vgh_sgs_flux = flip(vgh_sgs_flux,3);
            ugh_sgs_flux = flip(ugh_sgs_flux,3);
                        
            sum_vghsgs = sum(vgh_sgs,3);
            sum_ughsgs = sum(ugh_sgs,3);              
            sum_hfxy = sum(h_fxy,3);
            for l=1:xsize_f
                for m=1:ysize_f
                    vghb_sgs(l,m,:) = vgh_sgs(l,m,:) - (h_fxy(l,m,:) * (sum_vghsgs(l,m)/sum_hfxy(l,m)));  %%% remove the barotropic part
                    ughb_sgs(l,m,:) = ugh_sgs(l,m,:) - (h_fxy(l,m,:) * (sum_ughsgs(l,m)/sum_hfxy(l,m)));         
                end    
            end
            
            vghb_sgs_flux(:,:,1)=vghb_sgs(:,:,6);
            ughb_sgs_flux(:,:,1)=ughb_sgs(:,:,6);
            for k=2:zsize
                vghb_sgs_flux(:,:,k) =  ( vghb_sgs(:,:,zsize+1-k) ) + vghb_sgs_flux(:,:,k-1);
                ughb_sgs_flux(:,:,k) =  ( ughb_sgs(:,:,zsize+1-k) ) + ughb_sgs_flux(:,:,k-1);    
            end
            vghb_sgs_flux = flip(vghb_sgs_flux,3);
            ughb_sgs_flux = flip(ughb_sgs_flux,3);
            
            vghbp_sgs_flux(:,:,1) = ( h_fxy(:,:,6) .* (sum_vghsgs(:,:) ./ sum_hfxy(:,:)) );
            for k =2:zsize
                vghbp_sgs_flux(:,:,k) = (h_fxy(:,:,zsize+1-k) .* (sum_vghsgs(:,:) ./ sum_hfxy(:,:))) + vghbp_sgs_flux(:,:,k-1);    
            end
            vghbp_sgs_flux = flip(vghbp_sgs_flux,3);
            
            vghbp_sgs_flux(isnan(vghbp_sgs_flux))=0;
            vghb_sgs_flux(isnan(vghb_sgs_flux))=0;
            
            [dhdx,dhdy]=grad(h_fxy,DX_fxy,DY_fxy);
            [dedx,dedy]=grad(e_fxy,DX_fxy,DY_fxy);
            [dvdx,dvdy]=grad(vg_fxy,DX_fxy,DY_fxy);
            
                                    
            for j=1:ysize_f
                dhdyMO(:,j,:) = dhdy(:,j,:) - ( beta(j) * (h_fxy(:,j,:) ./ CORCOEF_f(j)) );
            end
            dhdyMO(isnan(dhdyMO))=0;

                        
            DEDX2(:,:,:) = dedx(:,:,:) .* dedx(:,:,:);
            DHDX2(:,:,:) = dhdx(:,:,:) .* dhdx(:,:,:);
            DEDY2(:,:,:) = dedy(:,:,:) .* dedy(:,:,:);
            DHDY2(:,:,:) = dhdy(:,:,:) .* dhdy(:,:,:);
            DHDY2MO(:,:,:) = dhdyMO(:,:,:) .* dhdyMO(:,:,:);
                        
            
            nom_sgs = ( (vgh_sgs_flux .* dedy) ) + ( (ugh_sgs_flux .* dedx) );
            nombS_sgs = ( (vghb_sgs_flux .* dedy) ) + ( (ughb_sgs_flux .* dedx) );
            denom_sgs = (DEDX2 + DEDY2);
            
            kappa_sgs = nom_sgs ./ denom_sgs;
                    
            nomh_sgs = ( (vgh_sgs .* dhdy) ./ ((DHDY2+DHDX2).^0.5) ) + ( (ugh_sgs .* dhdx) ./ ((DHDX2+DHDY2).^0.5) );
            denomh_sgs = (DHDX2 + DHDY2).^0.5; 
            
            nomh_sgsT1 = ( (vgh_sgs .* dhdy) ./ ((DHDY2).^0.5) );     
            denomh_sgsT1 = (DHDY2).^0.5; 
                        
            nomh_sgsT2 = ( (vgh_sgs .* dhdy) + (ugh_sgs .* dhdx));
            denomh_sgsT2 = (DHDX2 + DHDY2);
            nom_MixL_sgsT2 = VEL_sgs .* ( (vgh_sgs .* dhdy) + (ugh_sgs .* dhdx));
            denom_MixL_sgsT2 = (VEL_sgs .* VEL_sgs) .* denomh_sgsT2;
            Kappa_ed = (nomh_sgsT2 ./ denomh_sgsT2);
            
                      
            nomh_sgsMO = ( (vgh_sgs .* dhdyMO) ./ ((DHDY2MO+DHDX2).^0.5) ) + ( (ugh_sgs .* dhdx) ./ ((DHDX2+DHDY2MO).^0.5) );
            denomh_sgsMO = (DHDX2 + DHDY2MO).^0.5;
            
                                    
            nomh_sgsMOT2 = ( (vgh_sgs .* dhdyMO) + (ugh_sgs .* dhdx));
            denomh_sgsMOT2 = (DHDX2 + DHDY2MO);
            nom_MixL_sgsMOT2 = VEL_sgs .* ( (vgh_sgs .* dhdyMO) + (ugh_sgs .* dhdx));
            denom_MixL_sgsMOT2 = (VEL_sgs .* VEL_sgs) .* denomh_sgsMOT2;
            Kappa_edMO = (nomh_sgsMOT2 ./ denomh_sgsMOT2);
                       
            dhdynorm = dhdy ./ ((DHDY2).^0.5);
            dhdyMOnorm = dhdyMO ./ ((DHDY2MO).^0.5);
            
            
            [d2hdx2,d2hdxdy]=grad(dhdx,DX_fxy,DY_fxy);
            [d2edx2,d2edxdy]=grad(dedx,DX_fxy,DY_fxy);
            [d2hdydx,d2hdy2]=grad(dhdy,DX_fxy,DY_fxy);
            [d2edydx,d2edy2]=grad(dedy,DX_fxy,DY_fxy);
            [d2hdydxMO,d2hdy2MO]=grad(dhdyMO,DX_fxy,DY_fxy);
            
            [ddx_ughsgs,ddy_ughsgs] = grad(ugh_sgs,DX_fxy,DY_fxy);
            [ddx_vghsgs,ddy_vghsgs] = grad(vgh_sgs,DX_fxy,DY_fxy);
            
                        
            [ddx_ughfxy,ddy_ughfxy] = grad(ugh_fxy,DX_fxy,DY_fxy);
            [ddx_vghfxy,ddy_vghfxy] = grad(vgh_fxy,DX_fxy,DY_fxy);

            lhdhdy = dhdy(:,:,:) ./ h_fxy(:,:,:);
            lhdhdyMO = dhdyMO(:,:,:) ./ h_fxy(:,:,:);
            lhdhdy(isnan(lhdhdy))=0;
            lhdhdyMO(isnan(lhdhdyMO))=0;
            
            sum_div = ddx_ughsgs + ddy_vghsgs;
            sum_diff_h = d2hdy2 + d2hdx2;
            sum_diff_hMO = d2hdy2MO + d2hdx2;
            sum_diff_e  = d2edy2 + d2edx2;
            
            
                       
            sum_av = sum_av +1;
            eeta(:,:,:) = eeta(:,:,:) + ee(:,:,:);
            eta_fxy(:,:,:) = eta_fxy(:,:,:) + e_fxy(:,:,:);
            hta_fxy(:,:,:) = hta_fxy(:,:,:) + h_fxy(:,:,:);
            vgh_fxyta(:,:,:) = vgh_fxyta(:,:,:) + vgh_fxy(:,:,:);
            ugh_fxyta(:,:,:) = ugh_fxyta(:,:,:) + ugh_fxy(:,:,:);
            vg_fxyh_fxyta(:,:,:) = vg_fxyh_fxyta(:,:,:) + (vg_fxyh_fxy(:,:,:));           
            ug_fxyh_fxyta(:,:,:) = ug_fxyh_fxyta(:,:,:) + (ug_fxyh_fxy(:,:,:));
            vgh_sgsta(:,:,:) = vgh_sgsta(:,:,:) + vgh_sgs(:,:,:);
            ugh_sgsta(:,:,:) = ugh_sgsta(:,:,:) + ugh_sgs(:,:,:);          
            vgh_sgs0ta(:,:,:) = vgh_sgs0ta(:,:,:) + vgh_sgs0(:,:,:);
            ugh_sgs0ta(:,:,:) = ugh_sgs0ta(:,:,:) + ugh_sgs0(:,:,:);            
            vgh_sgs2ta(:,:,:) = vgh_sgs2ta(:,:,:) + vgh_sgs2(:,:,:);
            ugh_sgs2ta(:,:,:) = ugh_sgs2ta(:,:,:) + ugh_sgs2(:,:,:);
            vgh_sgsGMta(:,:,:) =  vgh_sgsGMta(:,:,:) + vgh_sgsGM(:,:,:);
            ugh_sgsGMta(:,:,:) =  ugh_sgsGMta(:,:,:) + ugh_sgsGM(:,:,:);
            vgh_sgs_fluxta(:,:,:) = vgh_sgs_fluxta(:,:,:) + vgh_sgs_flux(:,:,:);
            ugh_sgs_fluxta(:,:,:) = ugh_sgs_fluxta(:,:,:) + ugh_sgs_flux(:,:,:);
            vghb_sgsta(:,:,:) = vghb_sgsta(:,:,:) + vghb_sgs(:,:,:);           
            ughb_sgsta(:,:,:) = ughb_sgsta(:,:,:) + ughb_sgs(:,:,:);
            vghb_sgs_fluxta(:,:,:) = vghb_sgs_fluxta(:,:,:) + vghb_sgs_flux(:,:,:);
            vghbp_sgs_fluxta(:,:,:) = vghbp_sgs_fluxta(:,:,:) + vghbp_sgs_flux(:,:,:);
            ughb_sgs_fluxta(:,:,:) = ughb_sgs_fluxta(:,:,:) + ughb_sgs_flux(:,:,:);
            ddx_ughsgsta(:,:,:) = ddx_ughsgsta(:,:,:) + ddx_ughsgs(:,:,:);
            ddy_ughsgsta(:,:,:) = ddy_ughsgsta(:,:,:) + ddy_ughsgs(:,:,:);
            ddx_vghsgsta(:,:,:) = ddx_vghsgsta(:,:,:) + ddx_vghsgs(:,:,:);
            ddy_vghsgsta(:,:,:) = ddy_vghsgsta(:,:,:) + ddy_vghsgs(:,:,:);
            nom_sgsta(:,:,:) =nom_sgsta(:,:,:)  + nom_sgs(:,:,:);
            denom_sgsta(:,:,:) = denom_sgsta(:,:,:)  + denom_sgs(:,:,:);
            kappa_sgsta(:,:,:) = kappa_sgsta(:,:,:) + kappa_sgs(:,:,:);
            nomh_sgsta(:,:,:) =nomh_sgsta(:,:,:)  + nomh_sgs(:,:,:);
            nomh_sgsT1ta(:,:,:) =nomh_sgsT1ta(:,:,:)  + nomh_sgsT1(:,:,:);
            nomh_sgsT2ta(:,:,:) =nomh_sgsT2ta(:,:,:)  + nomh_sgsT2(:,:,:);
            nomh_sgsMOta(:,:,:) =nomh_sgsMOta(:,:,:)  + nomh_sgsMO(:,:,:);
            nomh_sgsMOT2ta(:,:,:) =nomh_sgsMOT2ta(:,:,:)  + nomh_sgsMOT2(:,:,:);
            nombS_sgsta(:,:,:) =nombS_sgsta(:,:,:)  + nombS_sgs(:,:,:);
            denomh_sgsta(:,:,:) = denomh_sgsta(:,:,:)  + denomh_sgs(:,:,:);
            denomh_sgsT1ta(:,:,:) = denomh_sgsT1ta(:,:,:)  + denomh_sgsT1(:,:,:);
            denomh_sgsT2ta(:,:,:) = denomh_sgsT2ta(:,:,:)  + denomh_sgsT2(:,:,:);
            denomh_sgsMOta(:,:,:) = denomh_sgsMOta(:,:,:)  + denomh_sgsMO(:,:,:);
            denomh_sgsMOT2ta(:,:,:) = denomh_sgsMOT2ta(:,:,:)  + denomh_sgsMOT2(:,:,:);
            Kappa_edta(:,:,:) = Kappa_edta(:,:,:)  + Kappa_ed(:,:,:);
            Kappa_edMOta(:,:,:) = Kappa_edMOta(:,:,:)  + Kappa_edMO(:,:,:);
            dhdyta(:,:,:) = dhdyta(:,:,:) + dhdy(:,:,:);
            dhdyMOta(:,:,:) = dhdyMOta(:,:,:) + dhdyMO(:,:,:); 
            lhdhdyta(:,:,:) = lhdhdyta(:,:,:) + lhdhdy(:,:,:);
            lhdhdyMOta(:,:,:) = lhdhdyMOta(:,:,:) + lhdhdyMO(:,:,:);
            dedyta(:,:,:) = dedyta(:,:,:) + dedy(:,:,:);
            dhdxta(:,:,:) = dhdxta(:,:,:) + dhdx(:,:,:);
            dedxta(:,:,:) = dedxta(:,:,:) + dedx(:,:,:);
            d2hdy2ta(:,:,:) = d2hdy2ta(:,:,:) + d2hdy2(:,:,:);
            d2hdy2MOta(:,:,:) = d2hdy2MOta(:,:,:) + d2hdy2MO(:,:,:);
            d2edy2ta(:,:,:) = d2edy2ta(:,:,:) + d2edy2(:,:,:);
            d2hdx2ta(:,:,:) = d2hdx2ta(:,:,:) + d2hdx2(:,:,:);
            d2edx2ta(:,:,:) = d2edx2ta(:,:,:) + d2edx2(:,:,:);
            DHDY2ta(:,:,:)  = DHDY2ta(:,:,:) + DHDY2(:,:,:);
            DHDY2MOta(:,:,:) = DHDY2MOta(:,:,:) + DHDY2MO(:,:,:);
            DHDX2ta(:,:,:) = DHDX2ta(:,:,:) + DHDX2(:,:,:);
            vgvg_sgsta(:,:,:) = vgvg_sgsta(:,:,:) + vgvg_sgs(:,:,:);
            ugug_sgsta(:,:,:) = ugug_sgsta(:,:,:) + ugug_sgs(:,:,:);
            sum_divta(:,:,:)  = sum_divta(:,:,:) + sum_div(:,:,:);
            sum_diff_hta(:,:,:) = sum_diff_hta(:,:,:) + sum_diff_h(:,:,:);
            sum_diff_hMOta(:,:,:) = sum_diff_hMOta(:,:,:) + sum_diff_hMO(:,:,:);
            sum_diff_eta(:,:,:)   = sum_diff_eta(:,:,:) + sum_diff_e(:,:,:);
            nom_MixL_sgsT2ta(:,:,:) = nom_MixL_sgsT2ta(:,:,:) + nom_MixL_sgsT2(:,:,:);
            nom_MixL_sgsMOT2ta(:,:,:) = nom_MixL_sgsMOT2ta(:,:,:) + nom_MixL_sgsMOT2(:,:,:);
            denom_MixL_sgsT2ta(:,:,:) = denom_MixL_sgsT2ta(:,:,:) + denom_MixL_sgsT2(:,:,:);
            denom_MixL_sgsMOT2ta(:,:,:) = denom_MixL_sgsMOT2ta(:,:,:) + denom_MixL_sgsMOT2(:,:,:);
            VEL_sgsta(:,:,:) = VEL_sgsta(:,:,:) + VEL_sgs(:,:,:);
            
        end
end



eta_fxy = eta_fxy / sum_av;
eeta = eeta / sum_av;

[ugh_sgsta_ST,vgh_sgsta_ST,hta_fxy_ST,ugh_fxyta_ST,ug_fxyh_fxyta_ST,vgh_fxyta_ST,vg_fxyh_fxyta_ST]=SGSflux_cheap_NV(eeta,CORCOEF,gprime,DX,DY,nfxy,eps,1); %%%%% crops 1 is the best so far
[ugug_sgsta_ST,vgvg_sgsta_ST]=SGSeddvis_cheap_NV(eeta,CORCOEF,gprime,DX,DY,nfxy,eps,1);

ugh_sgsta_ST(isnan(ugh_sgsta_ST))=0;
vgh_sgsta_ST(isnan(vgh_sgsta_ST))=0;
hta_fxy_ST(isnan(hta_fxy_ST))=0;
ugh_fxyta_ST(isnan(ugh_fxyta_ST))=0;
ug_fxyh_fxyta_ST(isnan(ug_fxyh_fxyta_ST))=0;
vgh_fxyta_ST(isnan(vgh_fxyta_ST))=0;
vg_fxyh_fxyta_ST(isnan(vg_fxyh_fxyta_ST))=0;

ugug_sgsta_ST(isnan(ugug_sgsta_ST))=0;
vgvg_sgsta_ST(isnan(vgvg_sgsta_ST))=0;


DHDY2ta = DHDY2ta / sum_av;
DHDY2MOta = DHDY2MOta / sum_av;
DHDX2ta = DHDX2ta / sum_av;

vta_fxy = vta_fxy / sum_av;
uta_fxy = uta_fxy / sum_av;


vgh_fxyta = vgh_fxyta / sum_av;
ugh_fxyta = ugh_fxyta / sum_av;
vg_fxyh_fxyta = vg_fxyh_fxyta / sum_av;
ug_fxyh_fxyta = ug_fxyh_fxyta / sum_av;
vgh_sgsta = vgh_sgsta / sum_av;
ugh_sgsta = ugh_sgsta / sum_av;
vgh_sgs0ta = vgh_sgs0ta / sum_av;
ugh_sgs0ta = ugh_sgs0ta / sum_av;
vgh_sgs2ta = vgh_sgs2ta / sum_av;
ugh_sgs2ta = ugh_sgs2ta / sum_av;
vgh_sgsGMta = vgh_sgsGMta / sum_av;
ugh_sgsGMta = ugh_sgsGMta / sum_av;
ugh_sgs_fluxta = ugh_sgs_fluxta /sum_av;
vgh_sgs_fluxta = vgh_sgs_fluxta /sum_av;
vghb_sgsta = vghb_sgsta / sum_av;
vghbp_sgs_fluxta = vghbp_sgs_fluxta / sum_av;
ughb_sgsta = ughb_sgsta / sum_av;
vghb_sgs_fluxta = vghb_sgs_fluxta /sum_av;
ughb_sgs_fluxta = ughb_sgs_fluxta /sum_av;
nom_sgsta = nom_sgsta/ sum_av;
denom_sgsta = denom_sgsta / sum_av;
kappa_sgsta = kappa_sgsta / sum_av;
nomh_sgsta = nomh_sgsta/ sum_av;
nomh_sgsT1ta = nomh_sgsT1ta/ sum_av;
nomh_sgsT2ta = nomh_sgsT2ta/ sum_av;
nomh_sgsMOta = nomh_sgsMOta/ sum_av;
nomh_sgsMOT2ta = nomh_sgsMOT2ta/ sum_av;
nombS_sgsta = nombS_sgsta/ sum_av;
denomh_sgsta = denomh_sgsta / sum_av;
denomh_sgsT1ta = denomh_sgsT1ta / sum_av;
denomh_sgsT2ta = denomh_sgsT2ta / sum_av;
denomh_sgsMOta = denomh_sgsMOta / sum_av;
denomh_sgsMOT2ta = denomh_sgsMOT2ta / sum_av;
Kappa_edta = Kappa_edta / sum_av;
Kappa_edMOta = Kappa_edMOta / sum_av;
dhdyta = dhdyta / sum_av;
dhdyMOta = dhdyMOta / sum_av;
lhdhdyta = lhdhdyta / sum_av;
lhdhdyMOta = lhdhdyMOta / sum_av;
dedyta = dedyta / sum_av;
dhdxta = dhdxta / sum_av;
dedxta = dedxta / sum_av;
d2hdy2ta = d2hdy2ta / sum_av;
d2hdy2MOta = d2hdy2MOta / sum_av;
d2edy2ta = d2edy2ta / sum_av;
d2hdx2ta = d2hdx2ta / sum_av;
d2edx2ta = d2edx2ta / sum_av;
vgvg_sgsta = vgvg_sgsta / sum_av;
ugug_sgsta = ugug_sgsta / sum_av;
ddx_ughsgsta = ddx_ughsgsta / sum_av;
ddx_vghsgsta = ddx_vghsgsta / sum_av;
ddy_ughsgsta = ddy_ughsgsta / sum_av;
ddy_vghsgsta = ddy_vghsgsta / sum_av;            
sum_divta = sum_divta / sum_av;
sum_diff_hta = sum_diff_hta / sum_av;
sum_diff_hMOta = sum_diff_hMOta/ sum_av;
sum_diff_eta = sum_diff_eta/ sum_av;
nom_MixL_sgsT2ta = nom_MixL_sgsT2ta / sum_av;
nom_MixL_sgsMOT2ta = nom_MixL_sgsMOT2ta / sum_av;
denom_MixL_sgsT2ta = denom_MixL_sgsT2ta / sum_av;
denom_MixL_sgsMOT2ta = denom_MixL_sgsMOT2ta / sum_av;
VEL_sgsta = VEL_sgsta / sum_av;
        
EKE_sgsta_ST = (ugug_sgsta_ST + vgvg_sgsta_ST);
VEL_sgsta_ST = real(EKE_sgsta_ST.^0.5);

dhdyMOta(isnan(dhdyMOta))=0;
dhdyta(isnan(dhdyta))=0;
lhdhdyMOta(isnan(lhdhdyMOta))=0;
lhdhdyta(isnan(lhdhdyta))=0;
eeta(isnan(eeta))=0;


vgh_sgsta_TR = zeros(xsize_f,ysize_f,zsize,'single');
ugh_sgsta_TR = zeros(xsize_f,ysize_f,zsize,'single');
ugug_sgsta_TR = zeros(xsize_f,ysize_f,zsize,'single');
vgvg_sgsta_TR = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsT2_TRta  = zeros(xsize_f,ysize_f,zsize,'single');
nomh_sgsMOT2_TRta  = zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsT2_TRta = zeros(xsize_f,ysize_f,zsize,'single');
nom_MixL_sgsMOT2_TRta = zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsT2_TRta = zeros(xsize_f,ysize_f,zsize,'single');
denom_MixL_sgsMOT2_TRta = zeros(xsize_f,ysize_f,zsize,'single');
VEL_sgsta_TR     = zeros(xsize_f,ysize_f,zsize,'single'); 

sum_av2 = 0;
for i= 1:tsize
        if (i>=TS)
            eee2 = ncread(openstr,'e',[1,1,1,i],[xsize,ysize,zsize+1,1]);

            eee2(eee2<-1.e33) = NaN;
            
            [ugh_sgs_TR,vgh_sgs_TR,h_fxy_TR,ugh_fxy_TR,ug_fxyh_fxy_TR,vgh_fxy_TR,vg_fxyh_fxy_TR]=SGSflux_cheap_NV(eee2,CORCOEF,gprime,DX,DY,nfxy,eps,1); %%% crops 1 is the best option 
            [ugug_sgs_TR,vgvg_sgs_TR]=SGSeddvis_cheap_NV(eee2,CORCOEF,gprime,DX,DY,nfxy,eps,1); 
            
            ugh_sgs_TR(isnan(ugh_sgs_TR))=0;
            vgh_sgs_TR(isnan(vgh_sgs_TR))=0;             
            h_fxy_TR(isnan(h_fxy_TR))=0;
            ugh_fxy_TR(isnan(ugh_fxy_TR))=0;
            ug_fxyh_fxy_TR(isnan(ug_fxyh_fxy_TR))=0;
            vgh_fxy_TR(isnan(vgh_fxy_TR))=0;
            ugug_sgs_TR(isnan(ugug_sgs_TR))=0;
            vgvg_sgs_TR(isnan(vgvg_sgs_TR))=0;
            
            ugh_sgs_TR = ugh_sgs_TR - ugh_sgsta_ST;
            vgh_sgs_TR = vgh_sgs_TR - vgh_sgsta_ST;
            ugug_sgs_TR = ugug_sgs_TR - ugug_sgsta_ST;
            vgvg_sgs_TR = vgvg_sgs_TR - vgvg_sgsta_ST;
            
            
            %EKE_sgs_TR = ( (ugug_sgs_TR .* ugug_sgs_TR) + (vgvg_sgs_TR .* vgvg_sgs_TR) ).^0.5;
            EKE_sgs_TR = ( ugug_sgs_TR  + vgvg_sgs_TR );
            VEL_sgs_TR = real(EKE_sgs_TR.^0.5);
            
            [dhdx,dhdy]=grad(h_fxy_TR,DX_fxy,DY_fxy);
                                    
            for j=1:ysize_f
                dhdyMO(:,j,:) = dhdy(:,j,:) - ( beta(j) * (h_fxy(:,j,:) ./ CORCOEF_f(j)) );
            end
            dhdyMO(isnan(dhdyMO))=0;

                                    
            nomh_sgsT2_TR = ( (vgh_sgs_TR .* dhdy) + (ugh_sgs_TR .* dhdx));
            nomh_sgsMOT2_TR = ( (vgh_sgs_TR .* dhdyMO) + (ugh_sgs_TR .* dhdx));
            nom_MixL_sgsT2_TR = VEL_sgs_TR .*  ((vgh_sgs_TR .* dhdy) + (ugh_sgs_TR .* dhdx));
            nom_MixL_sgsMOT2_TR = VEL_sgs_TR .*  ((vgh_sgs_TR .* dhdyMO) + (ugh_sgs_TR .* dhdx));
            denom_MixL_sgsT2_TR = (VEL_sgs_TR .* VEL_sgs_TR) .* ( (dhdx .* dhdx) + (dhdy .* dhdy) );
            denom_MixL_sgsMOT2_TR = (VEL_sgs_TR .* VEL_sgs_TR) .* ( (dhdx .* dhdx) + (dhdyMO .* dhdyMO) );
            
            sum_av2 = sum_av2 +1;
            vgh_sgsta_TR(:,:,:) = vgh_sgsta_TR(:,:,:) + vgh_sgs_TR(:,:,:);   
            ugh_sgsta_TR(:,:,:) = ugh_sgsta_TR(:,:,:) + ugh_sgs_TR(:,:,:);
            ugug_sgsta_TR(:,:,:) = ugug_sgsta_TR(:,:,:) + ugug_sgs_TR(:,:,:);
            vgvg_sgsta_TR(:,:,:) = vgvg_sgsta_TR(:,:,:) + vgvg_sgs_TR(:,:,:);
            nomh_sgsT2_TRta(:,:,:) = nomh_sgsT2_TRta(:,:,:) + nomh_sgsT2_TR(:,:,:);
            nomh_sgsMOT2_TRta(:,:,:) = nomh_sgsMOT2_TRta(:,:,:) + nomh_sgsMOT2_TR(:,:,:);
            nom_MixL_sgsT2_TRta(:,:,:) =  nom_MixL_sgsT2_TRta(:,:,:) + nom_MixL_sgsT2_TR(:,:,:);
            nom_MixL_sgsMOT2_TRta(:,:,:) = nom_MixL_sgsMOT2_TRta(:,:,:) + nom_MixL_sgsMOT2_TR(:,:,:);
            denom_MixL_sgsT2_TRta(:,:,:) = denom_MixL_sgsT2_TRta(:,:,:) + denom_MixL_sgsT2_TR(:,:,:);
            denom_MixL_sgsMOT2_TRta(:,:,:) = denom_MixL_sgsMOT2_TRta(:,:,:) + denom_MixL_sgsMOT2_TR(:,:,:);
            VEL_sgsta_TR(:,:,:) = VEL_sgsta_TR(:,:,:) + VEL_sgs_TR(:,:,:);

          
            
        end
end


nomh_sgsT2_TRta = nomh_sgsT2_TRta / sum_av2;
nomh_sgsMOT2_TRta = nomh_sgsMOT2_TRta / sum_av2;
nom_MixL_sgsT2_TRta = nom_MixL_sgsT2_TRta /  sum_av2;
nom_MixL_sgsMOT2_TRta = nom_MixL_sgsMOT2_TRta / sum_av2;
denom_MixL_sgsT2_TRta = denom_MixL_sgsT2_TRta / sum_av2;
denom_MixL_sgsMOT2_TRta = denom_MixL_sgsMOT2_TRta / sum_av2;
VEL_sgsta_TR = VEL_sgsta_TR / sum_av2;

vgh_sgsta_TR = vgh_sgsta - vgh_sgsta_ST; 
ugh_sgsta_TR = ugh_sgsta - ugh_sgsta_ST;


ugug_sgsta_TR = ugug_sgsta - ugug_sgsta_ST; 
vgvg_sgsta_TR = vgvg_sgsta - vgvg_sgsta_ST;




nom_MixL_sgsT2ta_ST = VEL_sgsta_ST .* ( (vgh_sgsta_ST .* dhdyta) + (ugh_sgsta_ST .* dhdxta));
nom_MixL_sgsMOT2ta_ST = VEL_sgsta_ST .* ( (vgh_sgsta_ST .* dhdyMOta) + (ugh_sgsta_ST .* dhdxta));
denom_MixL_sgsT2ta_ST =  (VEL_sgsta_ST .* VEL_sgsta_ST) .* (denomh_sgsT2ta);  
denom_MixL_sgsMOT2ta_ST =  (VEL_sgsta_ST .* VEL_sgsta_ST) .* (denomh_sgsMOT2ta);  

nomh_sgsta_ST = ( (vgh_sgsta_ST .* dhdyta) ./ ((DHDY2ta+DHDX2ta).^0.5) ) + ( (ugh_sgsta_ST .* dhdxta) ./ ((DHDX2ta+DHDY2ta).^0.5) );
nomh_sgsta_TR = ( (vgh_sgsta_TR .* dhdyta) ./ ((DHDY2ta+DHDX2ta).^0.5) ) + ( (ugh_sgsta_TR .* dhdxta) ./ ((DHDX2ta+DHDY2ta).^0.5) );
nomh_sgsT2ta_ST = ( (vgh_sgsta_ST .* dhdyta) + (ugh_sgsta_ST .* dhdxta));
nomh_sgsT2ta_TR = ( (vgh_sgsta_TR .* dhdyta) + (ugh_sgsta_TR .* dhdxta));
 
                      

nomh_sgsMOta_ST = ( (vgh_sgsta_ST .* dhdyMOta) ./ ((DHDY2MOta+DHDX2ta).^0.5) ) + ( (ugh_sgsta_ST .* dhdxta) ./ ((DHDX2ta+DHDY2MOta).^0.5) );
nomh_sgsMOta_TR = ( (vgh_sgsta_TR .* dhdyMOta) ./ ((DHDY2MOta+DHDX2ta).^0.5) ) + ( (ugh_sgsta_TR .* dhdxta) ./ ((DHDX2ta+DHDY2MOta).^0.5) );
nomh_sgsMOT2ta_ST = ( (vgh_sgsta_ST .* dhdyMOta) + (ugh_sgsta_ST .* dhdxta));
nomh_sgsMOT2ta_TR = ( (vgh_sgsta_TR .* dhdyMOta) + (ugh_sgsta_TR .* dhdxta));

MASKPV = ones(xsize,ysize,zsize,'single');

for l=1:zsize
    for n=1:ysize
        for m=1:xsize
            if (hta_fxy(m,n,l) < 10)
                MASKPV(m,n,l) =0;
            end
        end
    end
end

kappa =- nom_sgsta ./ denom_sgsta;

Y_f = y_f / (111.32 * 1000); %[deg]
X_f = x_f / (111.32 * 1000);  %[deg]

[ddx_ughsgsta_ST,ddy_ughsgsta_ST] = grad(ugh_sgsta_ST,DX_fxy,DY_fxy);
[ddx_ughsgsta_TR,ddy_ughsgsta_TR] = grad(ugh_sgsta_TR,DX_fxy,DY_fxy);
[ddx_ughsgsta_fluxta,ddy_ughsgsta_fluxta] = grad(ugh_sgs_fluxta,DX_fxy,DY_fxy);
[ddx_ughbsgs_fluxta,ddy_ughbsgs_fluxta] = grad(ughb_sgs_fluxta,DX_fxy,DY_fxy);

[ddx_vghsgsta_ST,ddy_vghsgsta_ST] = grad(vgh_sgsta_ST,DX_fxy,DY_fxy);
[ddx_vghsgsta_TR,ddy_vghsgsta_TR] = grad(vgh_sgsta_TR,DX_fxy,DY_fxy);
[ddx_vghsgsta_fluxta,ddy_vghsgsta_fluxta] = grad(vgh_sgs_fluxta,DX_fxy,DY_fxy);
[ddx_vghbsgs_fluxta,ddy_vghbsgs_fluxta] = grad(vghb_sgs_fluxta,DX_fxy,DY_fxy);
[ddx_kappa,ddy_kappa] = grad(kappa,DX_fxy,DY_fxy);


sum_div_ST = ddx_ughsgsta_ST + ddy_vghsgsta_ST;
sum_div_TR = ddx_ughsgsta_TR + ddy_vghsgsta_TR;
sum_flux_div = ddx_ughsgs_fluxta + ddy_vghsgs_fluxta;
sum_flux_divb = ddx_ughbsgs_fluxta + ddy_vghbsgs_fluxta;


toc

tic

eta_fxy_za = zeros(ysize_f,zsize,'single');
ugta_fxy_za = zeros(ysize_f,zsize,'single');
vgta_fxy_za = zeros(ysize_f,zsize,'single');
uta_fxy_za = zeros(ysize_f,zsize,'single');
vta_fxy_za = zeros(ysize_f,zsize,'single');
hta_fxy_za = zeros(ysize_f,zsize,'single');
DY_za = zeros(ysize_f,1,'single');
DX_za = zeros(xsize_f,1,'single');
vg_fxyh_fxyta_za = zeros(ysize_f,zsize,'single');
vgh_fxyta_za = zeros(ysize_f,zsize,'single');
vgh_sgsta_za = zeros(ysize_f,zsize,'single');
vgh_sgstaTR_za = zeros(ysize_f,zsize,'single');
vgh_sgstaST_za = zeros(ysize_f,zsize,'single');
ugh_fxyta_za= zeros(ysize_f,zsize,'single');
ug_fxyh_fxyta_za= zeros(ysize_f,zsize,'single');
ugh_sgsta_za = zeros(ysize_f,zsize,'single');
ugh_sgstaTR_za = zeros(ysize_f,zsize,'single');
ugh_sgstaST_za = zeros(ysize_f,zsize,'single');
vgh_sgs0ta_za = zeros(ysize_f,zsize,'single');
vgh_sgs2ta_za = zeros(ysize_f,zsize,'single');
vgh_sgsGMta_za = zeros(ysize_f,zsize,'single');
vghb_sgsta_za = zeros(ysize_f,zsize,'single');
vghb_sgs_fluxta_za = zeros(ysize_f,zsize,'single');
vghbp_sgs_fluxta_za = zeros(ysize_f,zsize,'single');
nom_sgsta_za = zeros(ysize_f,zsize,'single');
denom_sgsta_za = zeros(ysize_f,zsize,'single');
nomh_sgsta_za = zeros(ysize_f,zsize,'single');
denomh_sgsta_za = zeros(ysize_f,zsize,'single');
nomh_sgsT2ta_za = zeros(ysize_f,zsize,'single');
nomh_sgsT2taTR_za = zeros(ysize_f,zsize,'single');
nomh_sgsT2taST_za = zeros(ysize_f,zsize,'single');
nomh_sgsMOT2ta_za = zeros(ysize_f,zsize,'single');
nomh_sgsMOT2taTR_za = zeros(ysize_f,zsize,'single');
nomh_sgsMOT2taST_za = zeros(ysize_f,zsize,'single');
denomh_sgsT2ta_za = zeros(ysize_f,zsize,'single');
denomh_sgsMOT2ta_za = zeros(ysize_f,zsize,'single');
Kappa_edta_za    = zeros(ysize_f,zsize,'single');
Kappa_edMOta_za    = zeros(ysize_f,zsize,'single');
kappa_sgsta_za = zeros(ysize_f,zsize,'single');
kappa_za = zeros(ysize_f,zsize,'single');
ddx_ughsgsta_za = zeros(ysize_f,zsize,'single');
ddy_vghsgsta_za = zeros(ysize_f,zsize,'single');
ddy_kappa_za = zeros(ysize_f,zsize,'single');
dedyta_za = zeros(ysize_f,zsize,'single');
dhdyta_za = zeros(ysize_f,zsize,'single');
dhdyMOta_za = zeros(ysize_f,zsize,'single');
lhdhdyta_za = zeros(ysize_f,zsize,'single');
lhdhdyMOta_za = zeros(ysize_f,zsize,'single');
vgvg_sgsta_za = zeros(ysize_f,zsize,'single');
ugug_sgsta_za = zeros(ysize_f,zsize,'single');
vgvg_sgstaST_za = zeros(ysize_f,zsize,'single');
ugug_sgstaST_za = zeros(ysize_f,zsize,'single');
vgvg_sgstaTR_za = zeros(ysize_f,zsize,'single');
ugug_sgstaTR_za = zeros(ysize_f,zsize,'single');
nomh_sgsT2TRta_za = zeros(ysize_f,zsize,'single');
nomh_sgsMOT2TRta_za = zeros(ysize_f,zsize,'single');
denom_MixL_sgsT2ta_za = zeros(ysize_f,zsize,'single');
denom_MixL_sgsT2TRta_za = zeros(ysize_f,zsize,'single');
denom_MixL_sgsT2taST_za = zeros(ysize_f,zsize,'single');
denom_MixL_sgsMOT2ta_za = zeros(ysize_f,zsize,'single');
denom_MixL_sgsMOT2TRta_za = zeros(ysize_f,zsize,'single');
denom_MixL_sgsMOT2taST_za = zeros(ysize_f,zsize,'single');
VEL_sgsta_za = zeros(ysize_f,zsize,'single');
VEL_sgsTRta_za = zeros(ysize_f,zsize,'single');
VEL_sgstaST_za = zeros(ysize_f,zsize,'single');


MASK =ones(xsize_f,ysize_f,zsize,'single'); 
MASK_b =ones(xsize_f,ysize_f,zsize+1,'single'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MASKING LAND IN NEVERLAND 
sum_MAS = 0;   
for k=1:zsize             
    for i=1:xsize_f
        for j=1:ysize_f
            if (Y_f(j) >=-49.5 && X_f(i) <7.5)
                MASK(i,j,k) = 0;
                sum_MAS = sum_MAS +1;
            end
            if (Y_f(j) >=-49.5 && X_f(i) >82.5)
                MASK(i,j,k) = 0;
                sum_MAS = sum_MAS +1;
            end
            if (Y_f(j)>=-38) &&  (37.5<= X_f(i) &&  X_f(i) <=52.5)
                MASK(i,j,k) = 0;
                sum_MAS = sum_MAS +1;
            end
            if ( -65<=Y_f(j) && Y_f(j) <=-49.5) && (X_f(i) <=18)  
               MASK(i,j,k) = 0;
                sum_MAS = sum_MAS +1;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MASKING LAND IN NEVERLAND
eta_fxy = eta_fxy .* MASK;
ugta_fxy = ugta_fxy .* MASK;
vgta_fxy = vgta_fxy .* MASK;
vta_fxy = vta_fxy .* MASK;
uta_fxy = uta_fxy .* MASK;
hta_fxy = hta_fxy .* MASK;
vgh_sgsta = vgh_sgsta .* MASK;
vgh_sgsta_ST = vgh_sgsta_ST .* MASK;
vgh_sgsta_TR = vgh_sgsta_TR .* MASK;
vgh_sgs0ta = vgh_sgs0ta .* MASK;
vgh_sgs2ta = vgh_sgs2ta .* MASK;
vgh_sgsGMta = vgh_sgsGMta .* MASK;
vghb_sgsta = vghb_sgsta .* MASK;
vghb_sgs_fluxta = vghb_sgs_fluxta .* MASK;
vghbp_sgs_fluxta = vghbp_sgs_fluxta.* MASK;
nom_sgsta = nom_sgsta .* MASK;
denom_sgsta = denom_sgsta .* MASK;
nomh_sgsta = nomh_sgsta .* MASK;
nomh_sgsT1ta = nomh_sgsT1ta .* MASK;
nomh_sgsT2ta = nomh_sgsT2ta .* MASK;
nomh_sgsT2ta_TR = nomh_sgsT2ta_TR .* MASK;
nomh_sgsT2_TRta = nomh_sgsT2_TRta .* MASK;
nomh_sgsT2ta_ST = nomh_sgsT2ta_ST .* MASK;
nomh_sgsMOT2ta = nomh_sgsMOT2ta .* MASK;
nomh_sgsMOT2ta_TR = nomh_sgsMOT2ta_TR .* MASK;
nomh_sgsMOT2_TRta = nomh_sgsMOT2_TRta .* MASK;
nomh_sgsMOT2ta_ST = nomh_sgsMOT2ta_ST .* MASK;
denomh_sgsta = denomh_sgsta .* MASK;
denomh_sgsT1ta = denomh_sgsT1ta .* MASK;
denomh_sgsT2ta = denomh_sgsT2ta .* MASK;
denomh_sgsMOT2ta = denomh_sgsMOT2ta .* MASK;
kappa_sgsta = kappa_sgsta .* MASK;
kappa = kappa .* MASK;
Kappa_edta = Kappa_edta .* MASK;
Kappa_edMOta = Kappa_edMOta .* MASK;
ddx_ughsgsta = ddx_ughsgsta .* MASK;
ddy_vghsgsta = ddy_vghsgsta .* MASK;
ddy_kappa = ddy_kappa .* MASK;
dedyta = dedyta .* MASK;
dhdyta = dhdyta .* MASK;
dhdyMOta = dhdyMOta .* MASK;
vgvg_sgsta = vgvg_sgsta .* MASK;
ugug_sgsta = ugug_sgsta .* MASK;
vgvg_sgsta_ST = vgvg_sgsta_ST .* MASK;
ugug_sgsta_ST = ugug_sgsta_ST .* MASK;
vgvg_sgsta_TR = vgvg_sgsta_TR .* MASK;
ugug_sgsta_TR = ugug_sgsta_TR .* MASK;
ugh_sgsta    = ugh_sgsta .* MASK;
ugh_sgsta_ST = ugh_sgsta_ST .* MASK;
ugh_sgsta_TR = ugh_sgsta_TR .* MASK;
ugh_fxyta    = ugh_fxyta .* MASK;
ug_fxyh_fxyta = ug_fxyh_fxyta .* MASK;
vgh_fxyta    = vgh_fxyta .* MASK;
vg_fxyh_fxyta = vg_fxyh_fxyta .* MASK;
denom_MixL_sgsT2ta = denom_MixL_sgsT2ta .*  MASK;
denom_MixL_sgsT2_TRta= denom_MixL_sgsT2_TRta .*   MASK;
denom_MixL_sgsT2ta_ST  = denom_MixL_sgsT2ta_ST  .*  MASK;
denom_MixL_sgsMOT2ta = denom_MixL_sgsMOT2ta .*  MASK;
denom_MixL_sgsMOT2_TRta= denom_MixL_sgsMOT2_TRta .*   MASK;
denom_MixL_sgsMOT2ta_ST  = denom_MixL_sgsMOT2ta_ST  .*  MASK;
VEL_sgsta_TR  = VEL_sgsta_TR .*  MASK;
VEL_sgsta_ST  = VEL_sgsta_ST .*  MASK;
VEL_sgsta     = VEL_sgsta .* MASK;
 
lhdhdyta = (lhdhdyta .* MASK) .* MASKPV;
lhdhdyMOta = (lhdhdyMOta .* MASK) .* MASKPV;


sum_ef =0;
for j=1:ysize_f
    for k=1:zsize
        for i=1:xsize_f
            if (MASK(i,j,k) ==1)
               sum_ef = sum_ef +1;
               eta_fxy_za(j,k) = eta_fxy_za(j,k) + (eta_fxy(i,j,k));
               hta_fxy_za(j,k) = hta_fxy_za(j,k) + (hta_fxy(i,j,k));
               uta_fxy_za(j,k) = uta_fxy_za(j,k) + (uta_fxy(i,j,k));
               vta_fxy_za(j,k) = vta_fxy_za(j,k) + (vta_fxy(i,j,k));
               ugta_fxy_za(j,k) = ugta_fxy_za(j,k) + (ugta_fxy(i,j,k));
               vgta_fxy_za(j,k) = vgta_fxy_za(j,k) + (vgta_fxy(i,j,k));
               vg_fxyh_fxyta_za(j,k)  = vg_fxyh_fxyta_za(j,k) + (vg_fxyh_fxyta(i,j,k) .* DX_fxy(i,j)) ;
               vgh_fxyta_za(j,k)  = vgh_fxyta_za(j,k) + (vgh_fxyta(i,j,k)  .*  DX_fxy(i,j));
               vgh_sgsta_za(j,k) = vgh_sgsta_za(j,k) + (vgh_sgsta(i,j,k)  .*  DX_fxy(i,j));
               vgh_sgstaTR_za(j,k) = vgh_sgstaTR_za(j,k) + (vgh_sgsta_TR(i,j,k)  .*  DX_fxy(i,j));
               vgh_sgstaST_za(j,k) = vgh_sgstaST_za(j,k) + (vgh_sgsta_ST(i,j,k)  .*  DX_fxy(i,j));
               ugh_fxyta_za(j,k) = ugh_fxyta_za(j,k) + (ugh_fxyta(i,j,k)  );
               ug_fxyh_fxyta_za(j,k) = ug_fxyh_fxyta_za(j,k) + (ug_fxyh_fxyta(i,j,k));
               ugh_sgsta_za(j,k) = ugh_sgsta_za(j,k) + (ugh_sgsta(i,j,k)  );
               ugh_sgstaTR_za(j,k) = ugh_sgstaTR_za(j,k) + (ugh_sgsta_TR(i,j,k)  );
               ugh_sgstaST_za(j,k) = ugh_sgstaST_za(j,k) + (ugh_sgsta_ST(i,j,k)  );
               vgh_sgs0ta_za(j,k) = vgh_sgs0ta_za(j,k) + (vgh_sgs0ta(i,j,k)  .*  DX_fxy(i,j));
               vgh_sgs2ta_za(j,k) = vgh_sgs2ta_za(j,k) + (vgh_sgs2ta(i,j,k)  .*  DX_fxy(i,j));
               vgh_sgsGMta_za(j,k) = vgh_sgsGMta_za(j,k) + (vgh_sgsGMta(i,j,k) .*  DX_fxy(i,j));
               vghb_sgsta_za(j,k) = vghb_sgsta_za(j,k) + (vghb_sgsta(i,j,k)  .*  DX_fxy(i,j));
               vghb_sgs_fluxta_za(j,k) = vghb_sgs_fluxta_za(j,k) + (vghb_sgs_fluxta(i,j,k)  .*  DX_fxy(i,j));
               vghbp_sgs_fluxta_za(j,k) = vghbp_sgs_fluxta_za(j,k) + (vghbp_sgs_fluxta(i,j,k)  .*  DX_fxy(i,j));
               nom_sgsta_za(j,k) = nom_sgsta_za(j,k) + (nom_sgsta(i,j,k)  .*  DX_fxy(i,j));
               denom_sgsta_za(j,k) = denom_sgsta_za(j,k) + (denom_sgsta(i,j,k) .* DX_fxy(i,j));
               nomh_sgsta_za(j,k) = nomh_sgsta_za(j,k) + (nomh_sgsta(i,j,k));
               denomh_sgsta_za(j,k) = denomh_sgsta_za(j,k) + (denomh_sgsta(i,j,k));               
               nomh_sgsT2ta_za(j,k) = nomh_sgsT2ta_za(j,k) + (nomh_sgsT2ta(i,j,k));
               nomh_sgsT2taTR_za(j,k) = nomh_sgsT2taTR_za(j,k) + (nomh_sgsT2ta_TR(i,j,k));
               nomh_sgsT2TRta_za(j,k) = nomh_sgsT2TRta_za(j,k) + (nomh_sgsT2_TRta(i,j,k));
               nomh_sgsT2taST_za(j,k) = nomh_sgsT2taST_za(j,k) + (nomh_sgsT2ta_ST(i,j,k));
               nomh_sgsMOT2ta_za(j,k) = nomh_sgsMOT2ta_za(j,k) + (nomh_sgsMOT2ta(i,j,k));
               nomh_sgsMOT2taTR_za(j,k) = nomh_sgsMOT2taTR_za(j,k) + (nomh_sgsMOT2ta_TR(i,j,k));
               nomh_sgsMOT2TRta_za(j,k) = nomh_sgsMOT2TRta_za(j,k) + (nomh_sgsMOT2_TRta(i,j,k));
               nomh_sgsMOT2taST_za(j,k) = nomh_sgsMOT2taST_za(j,k) + (nomh_sgsMOT2ta_ST(i,j,k));
               denomh_sgsT2ta_za(j,k) = denomh_sgsT2ta_za(j,k) + (denomh_sgsT2ta(i,j,k));
               denomh_sgsMOT2ta_za(j,k) = denomh_sgsMOT2ta_za(j,k) + (denomh_sgsMOT2ta(i,j,k));
               Kappa_edta_za(j,k) = Kappa_edta_za(j,k) + (Kappa_edta(i,j,k));
               Kappa_edMOta_za(j,k) = Kappa_edMOta_za(j,k) + (Kappa_edMOta(i,j,k));
               kappa_sgsta_za(j,k) = kappa_sgsta_za(j,k) + (kappa_sgsta(i,j,k)  .*  DX_fxy(i,j));
               kappa_za(j,k) = kappa_za(j,k) + (kappa(i,j,k)  .*  DX_fxy(i,j));
               ddx_ughsgsta_za(j,k) = ddx_ughsgsta_za(j,k) + (ddx_ughsgsta(i,j,k)  .*  DX_fxy(i,j));
               ddy_vghsgsta_za(j,k) = ddy_vghsgsta_za(j,k) + (ddy_vghsgsta(i,j,k)  .*  DX_fxy(i,j));
               ddy_kappa_za(j,k) = ddy_kappa_za(j,k) + (ddy_kappa(i,j,k)  .*  DX_fxy(i,j));
               dedyta_za(j,k)   = dedyta_za(j,k)  + ( dedyta(i,j,k));
               dhdyta_za(j,k)   = dhdyta_za(j,k)  + ( dhdyta(i,j,k));
               dhdyMOta_za(j,k)   = dhdyMOta_za(j,k)  + ( dhdyMOta(i,j,k));
               lhdhdyta_za(j,k)   = lhdhdyta_za(j,k)  + ( lhdhdyta(i,j,k));
               lhdhdyMOta_za(j,k)   = lhdhdyMOta_za(j,k)  + ( lhdhdyMOta(i,j,k));
               DY_za(j)        = DY_za(j) + ( DY_fxy(i,j));
               DX_za(j)        = DX_za(j) + ( DX_fxy(i,j));
               vgvg_sgsta_za(j,k) = vgvg_sgsta_za(j,k) + (vgvg_sgsta(i,j,k).*hta_fxy(i,j,k));
               ugug_sgsta_za(j,k) = ugug_sgsta_za(j,k) + (ugug_sgsta(i,j,k).*hta_fxy(i,j,k));
               vgvg_sgstaST_za(j,k) = vgvg_sgstaST_za(j,k) + (vgvg_sgsta_ST(i,j,k).*hta_fxy(i,j,k));
               ugug_sgstaST_za(j,k) = ugug_sgstaST_za(j,k) + (ugug_sgsta_ST(i,j,k).*hta_fxy(i,j,k));
               vgvg_sgstaTR_za(j,k) = vgvg_sgstaTR_za(j,k) + (vgvg_sgsta_TR(i,j,k).*hta_fxy(i,j,k));
               ugug_sgstaTR_za(j,k) = ugug_sgstaTR_za(j,k) + (ugug_sgsta_TR(i,j,k).*hta_fxy(i,j,k));
               denom_MixL_sgsT2ta_za(j,k)  = denom_MixL_sgsT2ta_za(j,k)  + denom_MixL_sgsT2ta(i,j,k);
               denom_MixL_sgsT2TRta_za(j,k) = denom_MixL_sgsT2TRta_za(j,k)     +  denom_MixL_sgsT2_TRta(i,j,k);
               denom_MixL_sgsT2taST_za(j,k) = denom_MixL_sgsT2taST_za(j,k) +  denom_MixL_sgsT2ta_ST(i,j,k);
               denom_MixL_sgsMOT2ta_za(j,k)  = denom_MixL_sgsMOT2ta_za(j,k)  + denom_MixL_sgsMOT2ta(i,j,k);
               denom_MixL_sgsMOT2TRta_za(j,k) = denom_MixL_sgsMOT2TRta_za(j,k)     +  denom_MixL_sgsMOT2_TRta(i,j,k);
               denom_MixL_sgsMOT2taST_za(j,k) = denom_MixL_sgsMOT2taST_za(j,k) +  denom_MixL_sgsMOT2ta_ST(i,j,k);
               VEL_sgsta_za(j,k)  = VEL_sgsta_za(j,k)  +  VEL_sgsta(i,j,k);
               VEL_sgsTRta_za(j,k) = VEL_sgsTRta_za(j,k) + VEL_sgsta_TR(i,j,k);
               VEL_sgstaST_za(j,k)  = VEL_sgstaST_za(j,k) + VEL_sgsta_ST(i,j,k);
            end
        end
    end
end

sum_ef = sum_ef /(ysize_f*(zsize));
sum_eft = sum_ef;
sum_ef = floor(sum_ef);

nomh_sgsta_za = nomh_sgsta_za /sum_ef;
denomh_sgsta_za = denomh_sgsta_za /sum_ef;
Kappa_c = nomh_sgsta_za ./ denomh_sgsta_za;
Kappa_c(isnan(Kappa_c)) = 0;
Kappa_c_VM = mean(Kappa_c(:,:),1);

nomh_sgsT2ta_za = nomh_sgsT2ta_za /sum_ef;
nomh_sgsT2taTR_za = nomh_sgsT2taTR_za /sum_ef;
nomh_sgsT2TRta_za = nomh_sgsT2TRta_za /sum_ef;
nomh_sgsT2taST_za = nomh_sgsT2taST_za /sum_ef;
nomh_sgsMOT2ta_za = nomh_sgsMOT2ta_za /sum_ef;
nomh_sgsMOT2taTR_za = nomh_sgsMOT2taTR_za /sum_ef;
nomh_sgsMOT2TRta_za = nomh_sgsMOT2TRta_za /sum_ef;
nomh_sgsMOT2taST_za = nomh_sgsMOT2taST_za /sum_ef;
denomh_sgsT2ta_za = denomh_sgsT2ta_za /sum_ef;
denomh_sgsMOT2ta_za = denomh_sgsMOT2ta_za /sum_ef;
denom_MixL_sgsT2ta_za = denom_MixL_sgsT2ta_za /sum_ef;
denom_MixL_sgsT2TRta_za = denom_MixL_sgsT2TRta_za /sum_ef; 
denom_MixL_sgsT2taST_za = denom_MixL_sgsT2taST_za /sum_ef;
denom_MixL_sgsMOT2ta_za = denom_MixL_sgsMOT2ta_za /sum_ef;
denom_MixL_sgsMOT2TRta_za = denom_MixL_sgsMOT2TRta_za /sum_ef; 
denom_MixL_sgsMOT2taST_za = denom_MixL_sgsMOT2taST_za /sum_ef;
VEL_sgsta_za = VEL_sgsta_za /sum_ef;
VEL_sgsTRta_za = VEL_sgsTRta_za /sum_ef;
VEL_sgstaST_za = VEL_sgstaST_za /sum_ef; 
Kappa_a = nomh_sgsT2ta_za ./ denomh_sgsT2ta_za;
Kappa_aTR = nomh_sgsT2taTR_za ./ denomh_sgsT2ta_za;
Kappa_aTRta = nomh_sgsT2TRta_za ./ denomh_sgsT2ta_za;
Kappa_aST = nomh_sgsT2taST_za ./ denomh_sgsT2ta_za;
KappaMO_a = nomh_sgsMOT2ta_za ./ denomh_sgsMOT2ta_za;
KappaMO_aTR = nomh_sgsMOT2taTR_za ./ denomh_sgsMOT2ta_za;
KappaMO_aTRta = nomh_sgsMOT2TRta_za ./ denomh_sgsMOT2ta_za;
KappaMO_aST = nomh_sgsMOT2taST_za ./ denomh_sgsMOT2ta_za;


MIX_LENGTH = nomh_sgsT2ta_za ./ denom_MixL_sgsT2ta_za;
MIX_LENGTH_TR = nomh_sgsT2TRta_za ./ denom_MixL_sgsT2TRta_za;
MIX_LENGTH_ST = nomh_sgsT2taST_za ./ denom_MixL_sgsT2taST_za;
ML_kappa_a      = Kappa_a ./ VEL_sgsta_za;
ML_kappa_aTRta      = Kappa_aTRta ./ VEL_sgsTRta_za;
ML_kappa_aST      = Kappa_aST ./ VEL_sgstaST_za; 
MIX_MOLENGTH = nomh_sgsMOT2ta_za ./ denom_MixL_sgsMOT2ta_za;
MIX_MOLENGTH_TR = nomh_sgsMOT2TRta_za ./ denom_MixL_sgsMOT2TRta_za;
MIX_MOLENGTH_ST = nomh_sgsMOT2taST_za ./ denom_MixL_sgsMOT2taST_za;
ML_kappaMO_a      = KappaMO_a ./ VEL_sgsta_za;
ML_kappaMO_aTRta      = KappaMO_aTRta ./ VEL_sgsTRta_za;
ML_kappaMO_aST      = KappaMO_aST ./ VEL_sgstaST_za; 




Kappa_a(isnan(Kappa_a)) = 0;
Kappa_a_VM = mean(Kappa_a(:,:),1);

Kappa_b = Kappa_edta_za /sum_ef;
KappaMO_b = Kappa_edMOta_za /sum_ef;
Kappa_b(isnan(Kappa_b)) = 0;
Kappa_b_VM = mean(Kappa_b(:,:),1);

dhdyta_za = dhdyta_za / sum_ef;
dhdyMOta_za = dhdyMOta_za / sum_ef;
lhdhdyta_za = lhdhdyta_za / sum_ef;
lhdhdyMOta_za = lhdhdyMOta_za / sum_ef;
dedyta_za = dedyta_za / sum_ef;
eta_fxy_za = eta_fxy_za / sum_ef;
hta_fxy_za = hta_fxy_za / sum_ef;
ugta_fxy_za = ugta_fxy_za / sum_ef;
vgta_fxy_za = vgta_fxy_za / sum_ef;
uta_fxy_za = uta_fxy_za / sum_ef;
vta_fxy_za = vta_fxy_za / sum_ef;
DY_za     = DY_za / sum_ef;
DX_za     = DX_za / sum_ef;
ugug_sgsta_za = ugug_sgsta_za / sum_ef;
vgvg_sgsta_za = vgvg_sgsta_za / sum_ef;
ugug_sgstaST_za = ugug_sgstaST_za / sum_ef;
vgvg_sgstaST_za = vgvg_sgstaST_za / sum_ef;
ugug_sgstaTR_za = ugug_sgstaTR_za / sum_ef;
vgvg_sgstaTR_za = vgvg_sgstaTR_za / sum_ef;
ugh_fxyta_za     = ugh_fxyta_za  / sum_ef;         
ug_fxyh_fxyta_za = ug_fxyh_fxyta_za / sum_ef;              
ugh_sgsta_za     = ugh_sgsta_za / sum_ef;          
ugh_sgstaTR_za    = ugh_sgstaTR_za / sum_ef;           
ugh_sgstaST_za   = ugh_sgstaST_za / sum_ef;

ENG_sgstaza = (1/2) * ( (ugug_sgsta_za .* ugug_sgsta_za) + (vgvg_sgsta_za .* vgvg_sgsta_za) ).^(0.5);
ENG_sgstaza_TR = (1/2) * ( (ugug_sgstaTR_za .* ugug_sgstaTR_za) + (vgvg_sgstaTR_za .* vgvg_sgstaTR_za) ).^(0.5); 
ENG_sgstaza_ST = (1/2) * ( (ugug_sgstaST_za .* ugug_sgstaST_za) + (vgvg_sgstaST_za .* vgvg_sgstaST_za) ).^(0.5); 

LR_sgstaza = ( (VEL_sgsta_za).^(1/2) ) ./ (beta .^(1/2)); 
LR_sgstaza_TR = ( (VEL_sgsTRta_za).^(1/2) ) ./ (beta .^(1/2)); 
LR_sgstaza_ST = ( (VEL_sgstaST_za).^(1/2) ) ./ (beta .^(1/2)); 

dhdyMOta_za(isnan(dhdyMOta_za))=0;
dhdyta_za(isnan(dhdyta_za))=0;
lhdhdyMOta_za(isnan(lhdhdyMOta_za))=0;
lhdhdyta_za(isnan(lhdhdyta_za))=0;


%%%%%%%%%%%%%%%%%% Max or Min baroclinic transport
MAX_brc = max(vghb_sgs_fluxta_za');
MIN_brc = min(vghb_sgs_fluxta_za');

for m=1:ysize_f
    if abs(MAX_brc(m)) < abs(MIN_brc(m))
       MAX_brc(m) = MIN_brc(m);
    end
end

MAX_brt = max(vghbp_sgs_fluxta_za');
MIN_brt = min(vghbp_sgs_fluxta_za');

for m=1:ysize_f
    if abs(MAX_brt(m)) < abs(MIN_brt(m))
        MAX_brt(m) = MIN_brt(m);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Global average for masking large gradients near outcrops  
MASK_GA =ones(xsize_f,ysize_f,zsize,'single'); 
sum_MASK = 0;   
for k=1:zsize             
    for i=1:xsize_f
        for j=1:ysize_f
            if (denomh_sgsT2ta(i,j,k) > 2.5e-6)  %%%%%  Channel configurations 2.5e-6       
                MASK_GA(i,j,k) = 0;
                sum_MASK = sum_MASK +1;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Global average for masking large gradients near outcrops  

for i=1:xsize_f
    for j=1:ysize_f
        nom_sgsta_G(i,j,:) = nom_sgsta(i,j,:) *  AR_fxy(i,j);
        denom_sgsta_G(i,j,:) = denom_sgsta(i,j,:) *  AR_fxy(i,j);
        nomh_sgsta_G(i,j,:) = nomh_sgsta(i,j,:) *  AR_fxy(i,j);
        nomh_sgsT1ta_G(i,j,:) = nomh_sgsT1ta(i,j,:) *  AR_fxy(i,j);
        nomh_sgsT2ta_G(i,j,:) = nomh_sgsT2ta(i,j,:) *  AR_fxy(i,j);
        nomh_sgsMOT2ta_G(i,j,:) = nomh_sgsMOT2ta(i,j,:) *  AR_fxy(i,j);
        nomh_sgsta_TRG(i,j,:) = nomh_sgsta_TR(i,j,:) *  AR_fxy(i,j);
        nomh_sgsta_STG(i,j,:) = nomh_sgsta_ST(i,j,:) *  AR_fxy(i,j);
        nomh_sgsT2ta_STG(i,j,:) = nomh_sgsT2ta_ST(i,j,:) *  AR_fxy(i,j);
        nomh_sgsT2ta_TRG(i,j,:) = nomh_sgsT2ta_TR(i,j,:) *  AR_fxy(i,j);
        nomh_sgsT2_TRtaG(i,j,:) = nomh_sgsT2_TRta(i,j,:) *  AR_fxy(i,j);
        nomh_sgsMOT2ta_STG(i,j,:) = nomh_sgsMOT2ta_ST(i,j,:) *  AR_fxy(i,j);
        nomh_sgsMOT2ta_TRG(i,j,:) = nomh_sgsMOT2ta_TR(i,j,:) *  AR_fxy(i,j);
        nomh_sgsMOT2_TRtaG(i,j,:) = nomh_sgsMOT2_TRta(i,j,:) *  AR_fxy(i,j);
        nomh_sgsMOta_G(i,j,:) = nomh_sgsMOta(i,j,:) *  AR_fxy(i,j);
        nomh_sgsMOta_TRG(i,j,:) = nomh_sgsMOta_TR(i,j,:) *  AR_fxy(i,j);
        nomh_sgsMOta_STG(i,j,:) = nomh_sgsMOta_ST(i,j,:) *  AR_fxy(i,j);
        nombS_sgsta_G(i,j,:) = nombS_sgsta(i,j,:) *  AR_fxy(i,j);
        denomh_sgsta_G(i,j,:) = denomh_sgsta(i,j,:) *  AR_fxy(i,j);
        denomh_sgsT1ta_G(i,j,:) = denomh_sgsT1ta(i,j,:) *  AR_fxy(i,j);
        denomh_sgsT2ta_G(i,j,:) = denomh_sgsT2ta(i,j,:) *  AR_fxy(i,j);
        denomh_sgsMOT2ta_G(i,j,:) = denomh_sgsMOT2ta(i,j,:) *  AR_fxy(i,j);
        denomh_sgsMOta_G(i,j,:) = denomh_sgsMOta(i,j,:) *  AR_fxy(i,j);
        Kappa_edta_G(i,j,:) = Kappa_edta(i,j,:) *  AR_fxy(i,j);
        Kappa_edMOta_G(i,j,:) = Kappa_edMOta(i,j,:) *  AR_fxy(i,j);
        nom_MixL_sgsT2ta_G(i,j,:) = nom_MixL_sgsT2ta(i,j,:) *  AR_fxy(i,j);
        nom_MixL_sgsT2_TRtaG(i,j,:)= nom_MixL_sgsT2_TRta(i,j,:) *  AR_fxy(i,j);
        nom_MixL_sgsT2ta_STG(i,j,:)  = nom_MixL_sgsT2ta_ST(i,j,:)  *  AR_fxy(i,j);          
        nom_MixL_sgsMOT2ta_G(i,j,:) = nom_MixL_sgsMOT2ta(i,j,:) *  AR_fxy(i,j);
        nom_MixL_sgsMOT2_TRtaG(i,j,:)= nom_MixL_sgsMOT2_TRta(i,j,:) *  AR_fxy(i,j);
        nom_MixL_sgsMOT2ta_STG(i,j,:)  = nom_MixL_sgsMOT2ta_ST(i,j,:)  *  AR_fxy(i,j);
        denom_MixL_sgsT2ta_G(i,j,:) = denom_MixL_sgsT2ta(i,j,:) *  AR_fxy(i,j);
        denom_MixL_sgsT2_TRtaG(i,j,:)= denom_MixL_sgsT2_TRta(i,j,:) *  AR_fxy(i,j);
        denom_MixL_sgsT2ta_STG(i,j,:)  = denom_MixL_sgsT2ta_ST(i,j,:)  *  AR_fxy(i,j);          
        denom_MixL_sgsMOT2ta_G(i,j,:) = denom_MixL_sgsMOT2ta(i,j,:) *  AR_fxy(i,j);
        denom_MixL_sgsMOT2_TRtaG(i,j,:)= denom_MixL_sgsMOT2_TRta(i,j,:) *  AR_fxy(i,j);
        denom_MixL_sgsMOT2ta_STG(i,j,:)  = denom_MixL_sgsMOT2ta_ST(i,j,:)  *  AR_fxy(i,j);
        VEL_sgsta_TRG(i,j,:)  = VEL_sgsta_TR(i,j,:) *  AR_fxy(i,j);
        VEL_sgsta_STG(i,j,:)  = VEL_sgsta_ST(i,j,:) *  AR_fxy(i,j);
        VEL_sgsta_G(i,j,:)     = VEL_sgsta(i,j,:) *  AR_fxy(i,j);
 
    end
end
nom_sgsta_G = nom_sgsta_G .* MASK_GA .* MASK;
denom_sgsta_G = denom_sgsta_G .* MASK_GA .* MASK;
nombS_sgsta_G = nombS_sgsta_G .* MASK_GA .* MASK;
nomh_sgsta_G = nomh_sgsta_G .* MASK_GA .* MASK;
denomh_sgsta_G = denomh_sgsta_G .* MASK_GA .* MASK;
nomh_sgsT1ta_G = nomh_sgsT1ta_G .* MASK_GA .* MASK;
denomh_sgsT1ta_G = denomh_sgsT1ta_G .* MASK_GA .* MASK;
nomh_sgsT2ta_G = nomh_sgsT2ta_G .* MASK_GA .* MASK;
nomh_sgsT2ta_TRG = nomh_sgsT2ta_TRG .* MASK_GA .* MASK;
nomh_sgsT2ta_STG = nomh_sgsT2ta_STG .* MASK_GA .* MASK;
nomh_sgsMOT2ta_G = nomh_sgsMOT2ta_G .* MASK_GA .* MASK;
nomh_sgsMOT2ta_TRG = nomh_sgsMOT2ta_TRG .* MASK_GA .* MASK;
nomh_sgsMOT2ta_STG = nomh_sgsMOT2ta_STG .* MASK_GA .* MASK;
denomh_sgsT2ta_G = denomh_sgsT2ta_G .* MASK_GA .* MASK;
denomh_sgsMOT2ta_G = denomh_sgsMOT2ta_G .* MASK_GA .* MASK;
Kappa_edta_G =  Kappa_edta_G .*  MASK_GA .* MASK;
Kappa_edMOta_G =  Kappa_edMOta_G .*  MASK_GA .* MASK;
nomh_sgsT2_TRtaG = nomh_sgsT2_TRtaG .*  MASK_GA .* MASK;
nomh_sgsMOT2_TRtaG = nomh_sgsMOT2_TRtaG .*  MASK_GA .* MASK;
nom_MixL_sgsT2ta_G    =  nom_MixL_sgsT2ta_G .*  MASK_GA .* MASK;  
nom_MixL_sgsT2_TRtaG =  nom_MixL_sgsT2_TRtaG .*  MASK_GA .* MASK;
nom_MixL_sgsT2ta_STG =  nom_MixL_sgsT2ta_STG .*  MASK_GA .* MASK;
nom_MixL_sgsMOT2ta_G    =  nom_MixL_sgsMOT2ta_G .*  MASK_GA .* MASK;  
nom_MixL_sgsMOT2_TRtaG =  nom_MixL_sgsMOT2_TRtaG .*  MASK_GA .* MASK;
nom_MixL_sgsMOT2ta_STG =  nom_MixL_sgsMOT2ta_STG .*  MASK_GA .* MASK;
denom_MixL_sgsT2ta_G    =  denom_MixL_sgsT2ta_G .*  MASK_GA .* MASK;  
denom_MixL_sgsT2_TRtaG =  denom_MixL_sgsT2_TRtaG .*  MASK_GA .* MASK;
denom_MixL_sgsT2ta_STG =  denom_MixL_sgsT2ta_STG .*  MASK_GA .* MASK;
denom_MixL_sgsMOT2ta_G    =  denom_MixL_sgsMOT2ta_G .*  MASK_GA .* MASK;  
denom_MixL_sgsMOT2_TRtaG =  denom_MixL_sgsMOT2_TRtaG .*  MASK_GA .* MASK;
denom_MixL_sgsMOT2ta_STG =  denom_MixL_sgsMOT2ta_STG .*  MASK_GA .* MASK;
VEL_sgsta_TRG = VEL_sgsta_TRG.*  MASK_GA .* MASK;
VEL_sgsta_STG = VEL_sgsta_STG .*  MASK_GA .* MASK;
VEL_sgsta_G = VEL_sgsta_G .*  MASK_GA .* MASK;

AR_fxy_GA =0;
sum_areaav = 0;
for i=1:xsize_f
    for j=1:ysize_f
        sum_areaav = sum_areaav + 1;
        AR_fxy_GA = AR_fxy_GA + AR_fxy(i,j);
    end
end
AR_fxy_GA = AR_fxy_GA / sum_areaav;

nom_sgsta_GA =0;
denom_sgsta_GA=0;
nomh_sgsta_GA =0;
nomh_sgsT1ta_GA =0;
nomh_sgsT2ta_GA =0;
nomh_sgsMOT2ta_GA =0;
nomh_sgsta_TRGA =0;
nomh_sgsta_STGA =0;
nomh_sgsT2ta_TRGA =0;
nomh_sgsT2_TRtaGA =0;
nomh_sgsT2ta_STGA =0;
nomh_sgsMOta_GA =0;
nomh_sgsMOta_TRGA =0;
nomh_sgsMOta_STGA =0;
nomh_sgsMOT2ta_TRGA =0;
nomh_sgsMOT2_TRtaGA =0;
nomh_sgsMOT2ta_STGA =0;
nombS_sgsta_GA =0;
denomh_sgsta_GA=0;
denomh_sgsT1ta_GA=0;
denomh_sgsT2ta_GA=0;
denomh_sgsMOT2ta_GA=0;
denomh_sgsMOta_GA=0;
Kappa_edta_GA = 0;
Kappa_edMOta_GA = 0;
nom_MixL_sgsT2ta_GA = 0;   
nom_MixL_sgsT2_TRtaGA = 0;
nom_MixL_sgsT2ta_STGA = 0; 
nom_MixL_sgsMOT2ta_GA = 0;   
nom_MixL_sgsMOT2_TRtaGA = 0;
nom_MixL_sgsMOT2ta_STGA = 0;
denom_MixL_sgsT2ta_GA = 0;   
denom_MixL_sgsT2_TRtaGA = 0;
denom_MixL_sgsT2ta_STGA = 0; 
denom_MixL_sgsMOT2ta_GA = 0;   
denom_MixL_sgsMOT2_TRtaGA = 0;
denom_MixL_sgsMOT2ta_STGA = 0; 
VEL_sgsta_TRGA = 0;
VEL_sgsta_STGA = 0; 
VEL_sgsta_GA = 0; 


nomh_sgsT2ta_VGA =zeros(zsize,1);
denomh_sgsT2ta_VGA=zeros(zsize,1);
nomh_sgsta_VGA =zeros(zsize,1);
denomh_sgsta_VGA=zeros(zsize,1);
Kappa_edta_VGA=zeros(zsize,1);

nom_sgsta_GA1 =0;
denom_sgsta_GA1=0;
nomh_sgsta_GA1 =0;
nomh_sgsta_TRGA1 =0;
nomh_sgsta_STGA1 =0;
nomh_sgsMOta_GA1 =0;
nomh_sgsMOta_TRGA1 =0;
nomh_sgsMOta_STGA1 =0;
nombS_sgsta_GA1 =0;
denomh_sgsta_GA1=0;
denomh_sgsMOta_GA1=0;
nomh_sgsT2ta_GA1 =0;
nomh_sgsT2ta_STGA1 =0;
nomh_sgsT2ta_TRGA1 =0;
nomh_sgsT2_TRtaGA1 =0;
nomh_sgsMOT2ta_GA1 =0;
nomh_sgsMOT2ta_STGA1 =0;
nomh_sgsMOT2ta_TRGA1 =0;
nomh_sgsMOT2_TRtaGA1 =0;
denomh_sgsT2ta_GA1 = 0;
denomh_sgsMOT2ta_GA1 = 0;
nom_MixL_sgsT2ta_GA1 = 0;   
nom_MixL_sgsT2_TRtaGA1 = 0;
nom_MixL_sgsT2ta_STGA1 = 0; 
nom_MixL_sgsMOT2ta_GA1 = 0;   
nom_MixL_sgsMOT2_TRtaGA1 = 0;
nom_MixL_sgsMOT2ta_STGA1 = 0;
denom_MixL_sgsT2ta_GA1 = 0;   
denom_MixL_sgsT2_TRtaGA1 = 0;
denom_MixL_sgsT2ta_STGA1 = 0; 
denom_MixL_sgsMOT2ta_GA1 = 0;   
denom_MixL_sgsMOT2_TRtaGA1 = 0;
denom_MixL_sgsMOT2ta_STGA1 = 0;
VEL_sgsta_TRGA1 = 0;
VEL_sgsta_STGA1 = 0; 
VEL_sgsta_GA1 = 0; 

nom_sgsta_GA2 =0;
denom_sgsta_GA2=0;
nomh_sgsta_GA2 =0;
nomh_sgsta_TRGA2 =0;
nomh_sgsta_STGA2 =0;
nomh_sgsMOta_GA2 =0;
nomh_sgsMOta_TRGA2 =0;
nomh_sgsMOta_STGA2 =0;
nombS_sgsta_GA2 =0;
denomh_sgsta_GA2=0;
denomh_sgsMOta_GA2=0;
nomh_sgsT2ta_GA2 =0;
nomh_sgsT2ta_STGA2 =0;
nomh_sgsT2ta_TRGA2 =0;
nomh_sgsT2_TRtaGA2 =0;
nomh_sgsMOT2ta_GA2 =0;
nomh_sgsMOT2ta_STGA2 =0;
nomh_sgsMOT2ta_TRGA2 =0;
nomh_sgsMOT2_TRtaGA2 =0;
denomh_sgsT2ta_GA2 = 0;
denomh_sgsMOT2ta_GA2 = 0;
nom_MixL_sgsT2ta_GA2 = 0;   
nom_MixL_sgsT2_TRtaGA2 = 0;
nom_MixL_sgsT2ta_STGA2 = 0; 
nom_MixL_sgsMOT2ta_GA2 = 0;   
nom_MixL_sgsMOT2_TRtaGA2 = 0;
nom_MixL_sgsMOT2ta_STGA2 = 0;
denom_MixL_sgsT2ta_GA2 = 0;   
denom_MixL_sgsT2_TRtaGA2 = 0;
denom_MixL_sgsT2ta_STGA2 = 0; 
denom_MixL_sgsMOT2ta_GA2 = 0;   
denom_MixL_sgsMOT2_TRtaGA2 = 0;
denom_MixL_sgsMOT2ta_STGA2 = 0;
VEL_sgsta_TRGA2 = 0;
VEL_sgsta_STGA2 = 0; 
VEL_sgsta_GA2 = 0; 

nom_sgsta_GA3 =0;
denom_sgsta_GA3=0;
nomh_sgsta_GA3 =0;
nomh_sgsta_TRGA3 =0;
nomh_sgsta_STGA3 =0;
nomh_sgsMOta_GA3 =0;
nomh_sgsMOta_TRGA3 =0;
nomh_sgsMOta_STGA3 =0;
nombS_sgsta_GA3 =0;
denomh_sgsta_GA3=0;
denomh_sgsMOta_GA3=0;
nomh_sgsT2ta_GA3 =0;
nomh_sgsT2ta_STGA3 =0;
nomh_sgsT2ta_TRGA3 =0;
nomh_sgsT2_TRtaGA3 =0;
nomh_sgsMOT2ta_GA3 =0;
nomh_sgsMOT2ta_STGA3 =0;
nomh_sgsMOT2ta_TRGA3 =0;
nomh_sgsMOT2_TRtaGA3 =0;
denomh_sgsT2ta_GA3 = 0;
denomh_sgsMOT2ta_GA3 = 0;
denom_MixL_sgsT2ta_GA3 = 0;   
denom_MixL_sgsT2_TRtaGA3 = 0;
denom_MixL_sgsT2ta_STGA3 = 0; 
denom_MixL_sgsMOT2ta_GA3 = 0;   
denom_MixL_sgsMOT2_TRtaGA3 = 0;
denom_MixL_sgsMOT2ta_STGA3 = 0;
VEL_sgsta_TRGA3 = 0;
VEL_sgsta_STGA3 = 0; 
VEL_sgsta_GA3 = 0; 

nom_sgsta_G(isnan(nom_sgsta_G)) = 0;
denom_sgsta_G(isnan(denom_sgsta_G)) = 0;
nomh_sgsta_G(isnan(nomh_sgsta_G)) = 0;
nomh_sgsT1ta_G(isnan(nomh_sgsT1ta_G)) = 0;
nomh_sgsT2ta_G(isnan(nomh_sgsT2ta_G)) = 0;
nomh_sgsT2ta_TRG(isnan(nomh_sgsT2ta_TRG)) = 0;
nomh_sgsT2_TRtaG(isnan(nomh_sgsT2_TRtaG)) = 0;
nomh_sgsT2ta_STG(isnan(nomh_sgsT2ta_STG)) = 0;
nomh_sgsMOT2ta_G(isnan(nomh_sgsMOT2ta_G)) = 0;
nomh_sgsMOT2ta_TRG(isnan(nomh_sgsMOT2ta_TRG)) = 0;
nomh_sgsMOT2_TRtaG(isnan(nomh_sgsMOT2_TRtaG)) = 0;
nomh_sgsMOT2ta_STG(isnan(nomh_sgsMOT2ta_STG)) = 0;
nomh_sgsta_TRG(isnan(nomh_sgsta_TRG)) = 0;
nomh_sgsta_STG(isnan(nomh_sgsta_STG)) = 0;
nomh_sgsMOta_G(isnan(nomh_sgsMOta_G)) = 0;
nomh_sgsMOta_TRG(isnan(nomh_sgsMOta_TRG)) = 0;
nomh_sgsMOta_STG(isnan(nomh_sgsMOta_STG)) = 0;
nombS_sgsta_G(isnan(nombS_sgsta_G)) = 0;
denomh_sgsta_G(isnan(denomh_sgsta_G)) = 0;
denomh_sgsT1ta_G(isnan(denomh_sgsT1ta_G)) = 0;
denomh_sgsT2ta_G(isnan(denomh_sgsT2ta_G)) = 0;
denomh_sgsMOT2ta_G(isnan(denomh_sgsMOT2ta_G)) = 0;
denomh_sgsMOta_G(isnan(denomh_sgsMOta_G)) = 0;
Kappa_edta_G(isnan(Kappa_edta_G)) = 0;
Kappa_edMOta_G(isnan(Kappa_edMOta_G)) = 0;
nom_MixL_sgsT2ta_G(isnan(nom_MixL_sgsT2ta_G)) = 0;
nom_MixL_sgsT2_TRtaG(isnan(nom_MixL_sgsT2_TRtaG)) = 0;
nom_MixL_sgsT2ta_STG(isnan(nom_MixL_sgsT2ta_STG)) = 0;
nom_MixL_sgsMOT2ta_G(isnan(nom_MixL_sgsMOT2ta_G)) = 0;
nom_MixL_sgsMOT2_TRtaG(isnan(nom_MixL_sgsMOT2_TRtaG)) = 0;
nom_MixL_sgsMOT2ta_STG(isnan(nom_MixL_sgsMOT2ta_STG)) = 0;
denom_MixL_sgsT2ta_G(isnan(denom_MixL_sgsT2ta_G)) = 0;
denom_MixL_sgsT2_TRtaG(isnan(denom_MixL_sgsT2_TRtaG)) = 0;
denom_MixL_sgsT2ta_STG(isnan(denom_MixL_sgsT2ta_STG)) = 0;
denom_MixL_sgsMOT2ta_G(isnan(denom_MixL_sgsMOT2ta_G)) = 0;
denom_MixL_sgsMOT2_TRtaG(isnan(denom_MixL_sgsMOT2_TRtaG)) = 0;
denom_MixL_sgsMOT2ta_STG(isnan(denom_MixL_sgsMOT2ta_STG)) = 0;
VEL_sgsta_TRG(isnan(VEL_sgsta_TRG)) = 0;
VEL_sgsta_STG(isnan(VEL_sgsta_STG)) = 0;
VEL_sgsta_G(isnan(VEL_sgsta_G)) = 0;

sum_glob = 0;
sum_glob1 = 0;
sum_glob2 = 0;
sum_glob3 = 0;
for k=1:zsize
    for i=1:(xsize_f)  %%%
        for j=1:(ysize_f) %%
                sum_glob = sum_glob+1;
                nom_sgsta_GA = nom_sgsta_GA + nom_sgsta_G(i,j,k);
                denom_sgsta_GA = denom_sgsta_GA + denom_sgsta_G(i,j,k);
                nomh_sgsta_GA = nomh_sgsta_GA + nomh_sgsta_G(i,j,k);
                nomh_sgsT1ta_GA = nomh_sgsT1ta_GA + nomh_sgsT1ta_G(i,j,k);
                nomh_sgsT2ta_GA = nomh_sgsT2ta_GA + nomh_sgsT2ta_G(i,j,k);
                nomh_sgsT2ta_TRGA = nomh_sgsT2ta_TRGA + nomh_sgsT2ta_TRG(i,j,k);
                nomh_sgsT2_TRtaGA = nomh_sgsT2_TRtaGA + nomh_sgsT2_TRtaG(i,j,k);
                nomh_sgsT2ta_STGA = nomh_sgsT2ta_STGA + nomh_sgsT2ta_STG(i,j,k);
                nomh_sgsMOT2ta_GA = nomh_sgsMOT2ta_GA + nomh_sgsMOT2ta_G(i,j,k);
                nomh_sgsMOT2ta_TRGA = nomh_sgsMOT2ta_TRGA + nomh_sgsMOT2ta_TRG(i,j,k);
                nomh_sgsMOT2_TRtaGA = nomh_sgsMOT2_TRtaGA + nomh_sgsMOT2_TRtaG(i,j,k);
                nomh_sgsMOT2ta_STGA = nomh_sgsMOT2ta_STGA + nomh_sgsMOT2ta_STG(i,j,k);
                nomh_sgsta_TRGA = nomh_sgsta_TRGA + nomh_sgsta_TRG(i,j,k);
                nomh_sgsta_STGA = nomh_sgsta_STGA + nomh_sgsta_STG(i,j,k);
                nomh_sgsMOta_GA = nomh_sgsMOta_GA + nomh_sgsMOta_G(i,j,k);
                nomh_sgsMOta_TRGA = nomh_sgsMOta_TRGA + nomh_sgsMOta_TRG(i,j,k);
                nomh_sgsMOta_STGA = nomh_sgsMOta_STGA + nomh_sgsMOta_STG(i,j,k);
                nombS_sgsta_GA = nombS_sgsta_GA + nombS_sgsta_G(i,j,k);
                denomh_sgsta_GA = denomh_sgsta_GA + denomh_sgsta_G(i,j,k);
                denomh_sgsT1ta_GA = denomh_sgsT1ta_GA + denomh_sgsT1ta_G(i,j,k);
                denomh_sgsT2ta_GA = denomh_sgsT2ta_GA + denomh_sgsT2ta_G(i,j,k);
                denomh_sgsMOT2ta_GA = denomh_sgsMOT2ta_GA + denomh_sgsMOT2ta_G(i,j,k);
                denomh_sgsMOta_GA = denomh_sgsMOta_GA + denomh_sgsMOta_G(i,j,k);
                Kappa_edta_GA     = Kappa_edta_GA     + Kappa_edta_G(i,j,k);
                Kappa_edMOta_GA     = Kappa_edMOta_GA     + Kappa_edMOta_G(i,j,k);
                nomh_sgsT2ta_VGA(k) = nomh_sgsT2ta_VGA(k) + nomh_sgsT2ta_G(i,j,k);
                denomh_sgsT2ta_VGA(k) = denomh_sgsT2ta_VGA(k) + denomh_sgsT2ta_G(i,j,k);
                nomh_sgsta_VGA(k) = nomh_sgsta_VGA(k) + nomh_sgsta_G(i,j,k);
                denomh_sgsta_VGA(k) = denomh_sgsta_VGA(k) + denomh_sgsta_G(i,j,k);
                Kappa_edta_VGA(k)   = Kappa_edta_VGA(k) + Kappa_edta_G(i,j,k);   
                nom_MixL_sgsT2ta_GA = nom_MixL_sgsT2ta_GA + nom_MixL_sgsT2ta_G(i,j,k);
                nom_MixL_sgsT2_TRtaGA = nom_MixL_sgsT2_TRtaGA + nom_MixL_sgsT2_TRtaG(i,j,k);
                nom_MixL_sgsT2ta_STGA = nom_MixL_sgsT2ta_STGA + nom_MixL_sgsT2ta_STG(i,j,k);
                nom_MixL_sgsMOT2ta_GA = nom_MixL_sgsMOT2ta_GA + nom_MixL_sgsMOT2ta_G(i,j,k);
                nom_MixL_sgsMOT2_TRtaGA = nom_MixL_sgsMOT2_TRtaGA + nom_MixL_sgsMOT2_TRtaG(i,j,k);
                nom_MixL_sgsMOT2ta_STGA = nom_MixL_sgsMOT2ta_STGA + nom_MixL_sgsMOT2ta_STG(i,j,k);
                denom_MixL_sgsT2ta_GA = denom_MixL_sgsT2ta_GA + denom_MixL_sgsT2ta_G(i,j,k);
                denom_MixL_sgsT2_TRtaGA = denom_MixL_sgsT2_TRtaGA + denom_MixL_sgsT2_TRtaG(i,j,k);
                denom_MixL_sgsT2ta_STGA = denom_MixL_sgsT2ta_STGA + denom_MixL_sgsT2ta_STG(i,j,k);
                denom_MixL_sgsMOT2ta_GA = denom_MixL_sgsMOT2ta_GA + denom_MixL_sgsMOT2ta_G(i,j,k);
                denom_MixL_sgsMOT2_TRtaGA = denom_MixL_sgsMOT2_TRtaGA + denom_MixL_sgsMOT2_TRtaG(i,j,k);
                denom_MixL_sgsMOT2ta_STGA = denom_MixL_sgsMOT2ta_STGA + denom_MixL_sgsMOT2ta_STG(i,j,k);
                VEL_sgsta_TRGA = VEL_sgsta_TRGA + VEL_sgsta_TRG(i,j,k);
                VEL_sgsta_STGA = VEL_sgsta_STGA + VEL_sgsta_STG(i,j,k);
                VEL_sgsta_GA = VEL_sgsta_GA + VEL_sgsta_G(i,j,k);
                
                if(Y_f(j) >= -33.3)    %%% North
                    
                sum_glob1 = sum_glob1+1;
                nom_sgsta_GA1 = nom_sgsta_GA1 + nom_sgsta_G(i,j,k);
                denom_sgsta_GA1 = denom_sgsta_GA1 + denom_sgsta_G(i,j,k);
                nomh_sgsta_GA1 = nomh_sgsta_GA1 + nomh_sgsta_G(i,j,k);
                nomh_sgsta_TRGA1 = nomh_sgsta_TRGA1 + nomh_sgsta_TRG(i,j,k);
                nomh_sgsta_STGA1 = nomh_sgsta_STGA1 + nomh_sgsta_STG(i,j,k);
                nomh_sgsMOta_GA1 = nomh_sgsMOta_GA1 + nomh_sgsMOta_G(i,j,k);
                nomh_sgsMOta_TRGA1 = nomh_sgsMOta_TRGA1 + nomh_sgsMOta_TRG(i,j,k);
                nomh_sgsMOta_STGA1 = nomh_sgsMOta_STGA1 + nomh_sgsMOta_STG(i,j,k);
                nombS_sgsta_GA1 = nombS_sgsta_GA1 + nombS_sgsta_G(i,j,k);
                denomh_sgsta_GA1 = denomh_sgsta_GA1 + denomh_sgsta_G(i,j,k);
                denomh_sgsMOta_GA1 = denomh_sgsMOta_GA1 + denomh_sgsMOta_G(i,j,k);
                nomh_sgsT2ta_GA1 = nomh_sgsT2ta_GA1 + nomh_sgsT2ta_G(i,j,k);
                nomh_sgsT2ta_TRGA1 = nomh_sgsT2ta_TRGA1 + nomh_sgsT2ta_TRG(i,j,k);
                nomh_sgsT2_TRtaGA1 = nomh_sgsT2_TRtaGA1 + nomh_sgsT2_TRtaG(i,j,k);
                nomh_sgsT2ta_STGA1 = nomh_sgsT2ta_STGA1 + nomh_sgsT2ta_STG(i,j,k);
                nomh_sgsMOT2ta_GA1 = nomh_sgsMOT2ta_GA1 + nomh_sgsMOT2ta_G(i,j,k);
                nomh_sgsMOT2ta_TRGA1 = nomh_sgsMOT2ta_TRGA1 + nomh_sgsMOT2ta_TRG(i,j,k);
                nomh_sgsMOT2_TRtaGA1 = nomh_sgsMOT2_TRtaGA1 + nomh_sgsMOT2_TRtaG(i,j,k);
                nomh_sgsMOT2ta_STGA1 = nomh_sgsMOT2ta_STGA1 + nomh_sgsMOT2ta_STG(i,j,k);
                denomh_sgsT2ta_GA1 = denomh_sgsT2ta_GA1 + denomh_sgsT2ta_G(i,j,k);
                denomh_sgsMOT2ta_GA1 = denomh_sgsMOT2ta_GA1 + denomh_sgsMOT2ta_G(i,j,k); 
                nom_MixL_sgsT2ta_GA1 = nom_MixL_sgsT2ta_GA1 + nom_MixL_sgsT2ta_G(i,j,k);
                nom_MixL_sgsT2_TRtaGA1 = nom_MixL_sgsT2_TRtaGA1 + nom_MixL_sgsT2_TRtaG(i,j,k);
                nom_MixL_sgsT2ta_STGA1 = nom_MixL_sgsT2ta_STGA1 + nom_MixL_sgsT2ta_STG(i,j,k);
                nom_MixL_sgsMOT2ta_GA1 = nom_MixL_sgsMOT2ta_GA1 + nom_MixL_sgsMOT2ta_G(i,j,k);
                nom_MixL_sgsMOT2_TRtaGA1 = nom_MixL_sgsMOT2_TRtaGA1 + nom_MixL_sgsMOT2_TRtaG(i,j,k);
                nom_MixL_sgsMOT2ta_STGA1 = nom_MixL_sgsMOT2ta_STGA1 + nom_MixL_sgsMOT2ta_STG(i,j,k);
                denom_MixL_sgsT2ta_GA1 = denom_MixL_sgsT2ta_GA1 + denom_MixL_sgsT2ta_G(i,j,k);
                denom_MixL_sgsT2_TRtaGA1 = denom_MixL_sgsT2_TRtaGA1 + denom_MixL_sgsT2_TRtaG(i,j,k);
                denom_MixL_sgsT2ta_STGA1 = denom_MixL_sgsT2ta_STGA1 + denom_MixL_sgsT2ta_STG(i,j,k);
                denom_MixL_sgsMOT2ta_GA1 = denom_MixL_sgsMOT2ta_GA1 + denom_MixL_sgsMOT2ta_G(i,j,k);
                denom_MixL_sgsMOT2_TRtaGA1 = denom_MixL_sgsMOT2_TRtaGA1 + denom_MixL_sgsMOT2_TRtaG(i,j,k);
                denom_MixL_sgsMOT2ta_STGA1 = denom_MixL_sgsMOT2ta_STGA1 + denom_MixL_sgsMOT2ta_STG(i,j,k);
                VEL_sgsta_TRGA1 = VEL_sgsta_TRGA1 + VEL_sgsta_TRG(i,j,k);
                VEL_sgsta_STGA1 = VEL_sgsta_STGA1 + VEL_sgsta_STG(i,j,k);
                VEL_sgsta_GA1 = VEL_sgsta_GA1 + VEL_sgsta_G(i,j,k);
                
                
                        
                elseif (Y_f(j) < -50.00)  %%% South
                    
                sum_glob2 = sum_glob2+1;
                nom_sgsta_GA2 = nom_sgsta_GA2 + nom_sgsta_G(i,j,k);
                denom_sgsta_GA2 = denom_sgsta_GA2 + denom_sgsta_G(i,j,k);
                nomh_sgsta_GA2 = nomh_sgsta_GA2 + nomh_sgsta_G(i,j,k);
                nomh_sgsta_TRGA2 = nomh_sgsta_TRGA2 + nomh_sgsta_TRG(i,j,k);
                nomh_sgsta_STGA2 = nomh_sgsta_STGA2 + nomh_sgsta_STG(i,j,k);
                nomh_sgsMOta_GA2 = nomh_sgsMOta_GA2 + nomh_sgsMOta_G(i,j,k);
                nomh_sgsMOta_TRGA2 = nomh_sgsMOta_TRGA2 + nomh_sgsMOta_TRG(i,j,k);
                nomh_sgsMOta_STGA2 = nomh_sgsMOta_STGA2 + nomh_sgsMOta_STG(i,j,k);
                nombS_sgsta_GA2 = nombS_sgsta_GA2 + nombS_sgsta_G(i,j,k);
                denomh_sgsta_GA2 = denomh_sgsta_GA2 + denomh_sgsta_G(i,j,k);
                denomh_sgsMOta_GA2 = denomh_sgsMOta_GA2 + denomh_sgsMOta_G(i,j,k);
                nomh_sgsT2ta_GA2 = nomh_sgsT2ta_GA2 + nomh_sgsT2ta_G(i,j,k);
                nomh_sgsT2ta_TRGA2 = nomh_sgsT2ta_TRGA2 + nomh_sgsT2ta_TRG(i,j,k);
                nomh_sgsT2_TRtaGA2 = nomh_sgsT2_TRtaGA2 + nomh_sgsT2_TRtaG(i,j,k);
                nomh_sgsT2ta_STGA2 = nomh_sgsT2ta_STGA2 + nomh_sgsT2ta_STG(i,j,k);
                nomh_sgsMOT2ta_GA2 = nomh_sgsMOT2ta_GA2 + nomh_sgsMOT2ta_G(i,j,k);
                nomh_sgsMOT2ta_TRGA2 = nomh_sgsMOT2ta_TRGA2 + nomh_sgsMOT2ta_TRG(i,j,k);
                nomh_sgsMOT2_TRtaGA2 = nomh_sgsMOT2_TRtaGA2 + nomh_sgsMOT2_TRtaG(i,j,k);
                nomh_sgsMOT2ta_STGA2 = nomh_sgsMOT2ta_STGA2 + nomh_sgsMOT2ta_STG(i,j,k);
                denomh_sgsT2ta_GA2 = denomh_sgsT2ta_GA2 + denomh_sgsT2ta_G(i,j,k);
                denomh_sgsMOT2ta_GA2 = denomh_sgsMOT2ta_GA2 + denomh_sgsMOT2ta_G(i,j,k); 
                nom_MixL_sgsT2ta_GA2 = nom_MixL_sgsT2ta_GA2 + nom_MixL_sgsT2ta_G(i,j,k);
                nom_MixL_sgsT2_TRtaGA2 = nom_MixL_sgsT2_TRtaGA2 + nom_MixL_sgsT2_TRtaG(i,j,k);
                nom_MixL_sgsT2ta_STGA2 = nom_MixL_sgsT2ta_STGA2 + nom_MixL_sgsT2ta_STG(i,j,k);
                nom_MixL_sgsMOT2ta_GA2 = nom_MixL_sgsMOT2ta_GA2 + nom_MixL_sgsMOT2ta_G(i,j,k);
                nom_MixL_sgsMOT2_TRtaGA2 = nom_MixL_sgsMOT2_TRtaGA2 + nom_MixL_sgsMOT2_TRtaG(i,j,k);
                nom_MixL_sgsMOT2ta_STGA2 = nom_MixL_sgsMOT2ta_STGA2 + nom_MixL_sgsMOT2ta_STG(i,j,k);
                denom_MixL_sgsT2ta_GA2 = denom_MixL_sgsT2ta_GA2 + denom_MixL_sgsT2ta_G(i,j,k);
                denom_MixL_sgsT2_TRtaGA2 = denom_MixL_sgsT2_TRtaGA2 + denom_MixL_sgsT2_TRtaG(i,j,k);
                denom_MixL_sgsT2ta_STGA2 = denom_MixL_sgsT2ta_STGA2 + denom_MixL_sgsT2ta_STG(i,j,k);
                denom_MixL_sgsMOT2ta_GA2 = denom_MixL_sgsMOT2ta_GA2 + denom_MixL_sgsMOT2ta_G(i,j,k);
                denom_MixL_sgsMOT2_TRtaGA2 = denom_MixL_sgsMOT2_TRtaGA2 + denom_MixL_sgsMOT2_TRtaG(i,j,k);
                denom_MixL_sgsMOT2ta_STGA2 = denom_MixL_sgsMOT2ta_STGA2 + denom_MixL_sgsMOT2ta_STG(i,j,k);
                VEL_sgsta_TRGA2 = VEL_sgsta_TRGA2 + VEL_sgsta_TRG(i,j,k);
                VEL_sgsta_STGA2 = VEL_sgsta_STGA2 + VEL_sgsta_STG(i,j,k);
                VEL_sgsta_GA2 = VEL_sgsta_GA2 + VEL_sgsta_G(i,j,k);
                
                elseif ( -33.3 > Y_f(j) && Y_f(j) >= -50.00)  %%Central
                
                sum_glob3 = sum_glob3+1;
                nom_sgsta_GA3 = nom_sgsta_GA3 + nom_sgsta_G(i,j,k);
                denom_sgsta_GA3 = denom_sgsta_GA3 + denom_sgsta_G(i,j,k);
                nomh_sgsta_GA3 = nomh_sgsta_GA3 + nomh_sgsta_G(i,j,k);
                nomh_sgsta_TRGA3 = nomh_sgsta_TRGA3 + nomh_sgsta_TRG(i,j,k);
                nomh_sgsta_STGA3 = nomh_sgsta_STGA3 + nomh_sgsta_STG(i,j,k);
                nomh_sgsMOta_GA3 = nomh_sgsMOta_GA3 + nomh_sgsMOta_G(i,j,k);
                nomh_sgsMOta_TRGA3 = nomh_sgsMOta_TRGA3 + nomh_sgsMOta_TRG(i,j,k);
                nomh_sgsMOta_STGA3 = nomh_sgsMOta_STGA3 + nomh_sgsMOta_STG(i,j,k);
                nombS_sgsta_GA3 = nombS_sgsta_GA3 + nombS_sgsta_G(i,j,k);
                denomh_sgsta_GA3 = denomh_sgsta_GA3 + denomh_sgsta_G(i,j,k);
                denomh_sgsMOta_GA3 = denomh_sgsMOta_GA3 + denomh_sgsMOta_G(i,j,k);
                nomh_sgsT2ta_GA3 = nomh_sgsT2ta_GA3 + nomh_sgsT2ta_G(i,j,k);
                nomh_sgsT2ta_TRGA3 = nomh_sgsT2ta_TRGA3 + nomh_sgsT2ta_TRG(i,j,k);
                nomh_sgsT2_TRtaGA3 = nomh_sgsT2_TRtaGA3 + nomh_sgsT2_TRtaG(i,j,k);
                nomh_sgsT2ta_STGA3 = nomh_sgsT2ta_STGA3 + nomh_sgsT2ta_STG(i,j,k);
                nomh_sgsMOT2ta_GA3 = nomh_sgsMOT2ta_GA3 + nomh_sgsMOT2ta_G(i,j,k);
                nomh_sgsMOT2ta_TRGA3 = nomh_sgsMOT2ta_TRGA3 + nomh_sgsMOT2ta_TRG(i,j,k);
                nomh_sgsMOT2_TRtaGA3 = nomh_sgsMOT2_TRtaGA3 + nomh_sgsMOT2_TRtaG(i,j,k);
                nomh_sgsMOT2ta_STGA3 = nomh_sgsMOT2ta_STGA3 + nomh_sgsMOT2ta_STG(i,j,k);
                denomh_sgsT2ta_GA3 = denomh_sgsT2ta_GA3 + denomh_sgsT2ta_G(i,j,k);
                denomh_sgsMOT2ta_GA3 = denomh_sgsMOT2ta_GA3 + denomh_sgsMOT2ta_G(i,j,k);              
                denom_MixL_sgsT2ta_GA3 = denom_MixL_sgsT2ta_GA3 + denom_MixL_sgsT2ta_G(i,j,k);
                denom_MixL_sgsT2_TRtaGA3 = denom_MixL_sgsT2_TRtaGA3 + denom_MixL_sgsT2_TRtaG(i,j,k);
                denom_MixL_sgsT2ta_STGA3 = denom_MixL_sgsT2ta_STGA3 + denom_MixL_sgsT2ta_STG(i,j,k);
                denom_MixL_sgsMOT2ta_GA3 = denom_MixL_sgsMOT2ta_GA3 + denom_MixL_sgsMOT2ta_G(i,j,k);
                denom_MixL_sgsMOT2_TRtaGA3 = denom_MixL_sgsMOT2_TRtaGA3 + denom_MixL_sgsMOT2_TRtaG(i,j,k);
                denom_MixL_sgsMOT2ta_STGA3 = denom_MixL_sgsMOT2ta_STGA3 + denom_MixL_sgsMOT2ta_STG(i,j,k);
                VEL_sgsta_TRGA3 = VEL_sgsta_TRGA3 + VEL_sgsta_TRG(i,j,k);
                VEL_sgsta_STGA3 = VEL_sgsta_STGA3 + VEL_sgsta_STG(i,j,k);
                VEL_sgsta_GA3 = VEL_sgsta_GA3 + VEL_sgsta_G(i,j,k);
                
                
                end
         
        end
    end
end

KHTH_Ge =  -nom_sgsta_GA / denom_sgsta_GA;
KHTH_GebS = -nombS_sgsta_GA / denom_sgsta_GA;

KHTH_Gh = -nomh_sgsta_GA / denomh_sgsta_GA;
KHTH_GhT1 = -nomh_sgsT1ta_GA / denomh_sgsT1ta_GA;
KHTH_GhT2 = -nomh_sgsT2ta_GA / denomh_sgsT2ta_GA;
KHTH_TRGhT2 = -nomh_sgsT2ta_TRGA / denomh_sgsT2ta_GA;
KHTH_TRtaGhT2 = -nomh_sgsT2_TRtaGA / denomh_sgsT2ta_GA;
KHTH_STGhT2 = -nomh_sgsT2ta_STGA / denomh_sgsT2ta_GA;
KHTH_TRGh = -nomh_sgsta_TRGA / denomh_sgsta_GA;
KHTH_STGh = -nomh_sgsta_STGA / denomh_sgsta_GA;
KHTH_GhMO = -nomh_sgsMOta_GA / denomh_sgsMOta_GA;
KHTH_GhMOT2 = -nomh_sgsMOT2ta_GA / denomh_sgsMOT2ta_GA;
KHTH_TRGhMOT2 = -nomh_sgsMOT2ta_TRGA / denomh_sgsMOT2ta_GA;
KHTH_TRtaGhMOT2 = -nomh_sgsMOT2_TRtaGA / denomh_sgsMOT2ta_GA;
KHTH_STGhMOT2 = -nomh_sgsMOT2ta_STGA / denomh_sgsMOT2ta_GA;
KHTH_TRGhMO = -nomh_sgsMOta_TRGA / denomh_sgsMOta_GA;
KHTH_STGhMO = -nomh_sgsMOta_STGA / denomh_sgsMOta_GA;

VEL_sgsta_GA = (VEL_sgsta_GA / sum_glob) / AR_fxy_GA;
VEL_sgsta_TRGA = (VEL_sgsta_TRGA / sum_glob) / AR_fxy_GA;
VEL_sgsta_STGA = (VEL_sgsta_STGA / sum_glob) / AR_fxy_GA;

Lmix_ka = KHTH_GhT2 / VEL_sgsta_GA;
Lmix_kaTR = KHTH_TRtaGhT2 / VEL_sgsta_TRGA;
Lmix_kaST = KHTH_STGhT2 / VEL_sgsta_STGA;

Lmix     = -nom_MixL_sgsT2ta_GA / denom_MixL_sgsT2ta_GA;
Lmix_TR     = -nom_MixL_sgsT2_TRtaGA / denom_MixL_sgsT2_TRtaGA;
Lmix_ST     = -nom_MixL_sgsT2ta_STGA  / denom_MixL_sgsT2ta_STGA;

Lmix_MOka = KHTH_GhMOT2 / VEL_sgsta_GA;
Lmix_MOkaTR = KHTH_TRtaGhMOT2 / VEL_sgsta_TRGA;
Lmix_MOkaST = KHTH_STGhMOT2 / VEL_sgsta_STGA;

Lmix_MO     = -nom_MixL_sgsMOT2ta_GA / denom_MixL_sgsMOT2ta_GA;
Lmix_MOTR     = -nom_MixL_sgsMOT2_TRtaGA / denom_MixL_sgsMOT2_TRtaGA;
Lmix_MOST     = -nom_MixL_sgsMOT2ta_STGA  / denom_MixL_sgsMOT2ta_STGA;




KHTH_VGh = - nomh_sgsta_VGA ./ denomh_sgsta_VGA;
KHTH_VGhT2 = - nomh_sgsT2ta_VGA ./ denomh_sgsT2ta_VGA;

Kappa_edta_GA = Kappa_edta_GA / sum_glob;
Kappa_edMOta_GA = Kappa_edMOta_GA / sum_glob;
Kappa_edta_VGA = Kappa_edta_VGA / sum_areaav;
KAPPA_GH = - Kappa_edta_GA / AR_fxy_GA;
KAPPA_GHMO = - Kappa_edMOta_GA / AR_fxy_GA;
KAPPA_VGH = - Kappa_edta_VGA ./ AR_fxy_GA;



%%% $$$$$ North $$$$$

KHTH_Ge1 =  -nom_sgsta_GA1 / denom_sgsta_GA1;
KHTH_GebS1 = -nombS_sgsta_GA1 / denom_sgsta_GA1;

KHTH_Gh1 = -nomh_sgsta_GA1 / denomh_sgsta_GA1;
KHTH_TRGh1 = -nomh_sgsta_TRGA1 / denomh_sgsta_GA1;
KHTH_STGh1 = -nomh_sgsta_STGA1 / denomh_sgsta_GA1;
KHTH_GhMO1 = -nomh_sgsMOta_GA1 / denomh_sgsMOta_GA1;
KHTH_TRGhMO1 = -nomh_sgsMOta_TRGA1 / denomh_sgsMOta_GA1;
KHTH_STGhMO1 = -nomh_sgsMOta_STGA1 / denomh_sgsMOta_GA1;

KHTH_GhT21 = -nomh_sgsT2ta_GA1 / denomh_sgsT2ta_GA1;
KHTH_TRGhT21 = -nomh_sgsT2ta_TRGA1 / denomh_sgsT2ta_GA1;
KHTH_TRtaGhT21 = -nomh_sgsT2_TRtaGA1 / denomh_sgsT2ta_GA1;
KHTH_STGhT21 = -nomh_sgsT2ta_STGA1 / denomh_sgsT2ta_GA1;
KHTH_GhMOT21 = -nomh_sgsMOT2ta_GA1 / denomh_sgsMOT2ta_GA1;
KHTH_TRGhMOT21 = -nomh_sgsMOT2ta_TRGA1 / denomh_sgsMOT2ta_GA1;
KHTH_TRtaGhMOT21 = -nomh_sgsMOT2_TRtaGA1 / denomh_sgsMOT2ta_GA1;
KHTH_STGhMOT21 = -nomh_sgsMOT2ta_STGA1 / denomh_sgsMOT2ta_GA1;

VEL_sgsta_GA1 = (VEL_sgsta_GA1 / sum_glob1) / AR_fxy_GA;
VEL_sgsta_TRGA1 = (VEL_sgsta_TRGA1 / sum_glob1) / AR_fxy_GA;
VEL_sgsta_STGA1 = (VEL_sgsta_STGA1 / sum_glob1) / AR_fxy_GA;

Lmix_ka1 = KHTH_GhT21 / VEL_sgsta_GA1;
Lmix_kaTR1 = KHTH_TRtaGhT21 / VEL_sgsta_TRGA1;
Lmix_kaST1 = KHTH_STGhT21 / VEL_sgsta_STGA1;
Lmix1     = -nom_MixL_sgsT2ta_GA1 / denom_MixL_sgsT2ta_GA1;
Lmix_TR1     = -nom_MixL_sgsT2_TRtaGA1 / denom_MixL_sgsT2_TRtaGA1;
Lmix_ST1     = -nom_MixL_sgsT2ta_STGA1  / denom_MixL_sgsT2ta_STGA1;

Lmix_MOka1 = KHTH_GhMOT21 / VEL_sgsta_GA1;
Lmix_MOkaTR1 = KHTH_TRtaGhMOT21 / VEL_sgsta_TRGA1;
Lmix_MOkaST1 = KHTH_STGhMOT21 / VEL_sgsta_STGA1;
Lmix_MO1     = -nom_MixL_sgsMOT2ta_GA1 / denom_MixL_sgsMOT2ta_GA1;
Lmix_MOTR1     = -nom_MixL_sgsMOT2_TRtaGA1 / denom_MixL_sgsMOT2_TRtaGA1;
Lmix_MOST1     = -nom_MixL_sgsMOT2ta_STGA1  / denom_MixL_sgsMOT2ta_STGA1;


%%% $$$$$ South $$$$$

KHTH_Ge2 =  -nom_sgsta_GA2 / denom_sgsta_GA2;
KHTH_GebS2 = -nombS_sgsta_GA2 / denom_sgsta_GA2;

KHTH_Gh2 = -nomh_sgsta_GA2 / denomh_sgsta_GA2;
KHTH_TRGh2 = -nomh_sgsta_TRGA2 / denomh_sgsta_GA2;
KHTH_STGh2 = -nomh_sgsta_STGA2 / denomh_sgsta_GA2;
KHTH_GhMO2 = -nomh_sgsMOta_GA2 / denomh_sgsMOta_GA2;
KHTH_TRGhMO2 = -nomh_sgsMOta_TRGA2 / denomh_sgsMOta_GA2;
KHTH_STGhMO2 = -nomh_sgsMOta_STGA2 / denomh_sgsMOta_GA2;

KHTH_GhT22 = -nomh_sgsT2ta_GA2 / denomh_sgsT2ta_GA2;
KHTH_TRGhT22 = -nomh_sgsT2ta_TRGA2 / denomh_sgsT2ta_GA2;
KHTH_TRtaGhT22 = -nomh_sgsT2_TRtaGA2 / denomh_sgsT2ta_GA2;
KHTH_STGhT22 = -nomh_sgsT2ta_STGA2 / denomh_sgsT2ta_GA2;
KHTH_GhMOT22 = -nomh_sgsMOT2ta_GA2 / denomh_sgsMOT2ta_GA2;
KHTH_TRGhMOT22 = -nomh_sgsMOT2ta_TRGA2 / denomh_sgsMOT2ta_GA2;
KHTH_TRtaGhMOT22 = -nomh_sgsMOT2_TRtaGA2 / denomh_sgsMOT2ta_GA2;
KHTH_STGhMOT22 = -nomh_sgsMOT2ta_STGA2 / denomh_sgsMOT2ta_GA2;

VEL_sgsta_GA2 = (VEL_sgsta_GA2 / sum_glob2) / AR_fxy_GA;
VEL_sgsta_TRGA2 = (VEL_sgsta_TRGA2 / sum_glob2) / AR_fxy_GA;
VEL_sgsta_STGA2 = (VEL_sgsta_STGA2 / sum_glob2) / AR_fxy_GA;

Lmix_ka2 = KHTH_GhT22 / VEL_sgsta_GA2;
Lmix_kaTR2 = KHTH_TRtaGhT22 / VEL_sgsta_TRGA2;
Lmix_kaST2 = KHTH_STGhT22 / VEL_sgsta_STGA2;
Lmix2     = -nom_MixL_sgsT2ta_GA2 / denom_MixL_sgsT2ta_GA2;
Lmix_TR2     = -nom_MixL_sgsT2_TRtaGA2 / denom_MixL_sgsT2_TRtaGA2;
Lmix_ST2     = -nom_MixL_sgsT2ta_STGA2  / denom_MixL_sgsT2ta_STGA2;

Lmix_MOka2 = KHTH_GhMOT22 / VEL_sgsta_GA2;
Lmix_MOkaTR2 = KHTH_TRtaGhMOT22 / VEL_sgsta_TRGA2;
Lmix_MOkaST2 = KHTH_STGhMOT22 / VEL_sgsta_STGA2;
Lmix_MO2     = -nom_MixL_sgsMOT2ta_GA2 / denom_MixL_sgsMOT2ta_GA2;
Lmix_MOTR2     = -nom_MixL_sgsMOT2_TRtaGA2 / denom_MixL_sgsMOT2_TRtaGA2;
Lmix_MOST2     = -nom_MixL_sgsMOT2ta_STGA2  / denom_MixL_sgsMOT2ta_STGA2;


%%% $$$$$ Central $$$$$

KHTH_Ge3 =  -nom_sgsta_GA3 / denom_sgsta_GA3;
KHTH_GebS3 = -nombS_sgsta_GA3 / denom_sgsta_GA3;

KHTH_Gh3 = -nomh_sgsta_GA3 / denomh_sgsta_GA3;
KHTH_TRGh3 = -nomh_sgsta_TRGA3 / denomh_sgsta_GA3;
KHTH_STGh3 = -nomh_sgsta_STGA3 / denomh_sgsta_GA3;
KHTH_GhMO3 = -nomh_sgsMOta_GA3 / denomh_sgsMOta_GA3;
KHTH_TRGhMO3 = -nomh_sgsMOta_TRGA3 / denomh_sgsMOta_GA3;
KHTH_STGhMO3 = -nomh_sgsMOta_STGA3 / denomh_sgsMOta_GA3;

KHTH_GhT23 = -nomh_sgsT2ta_GA3 / denomh_sgsT2ta_GA3;
KHTH_TRGhT23 = -nomh_sgsT2ta_TRGA3 / denomh_sgsT2ta_GA3;
KHTH_TRtaGhT23 = -nomh_sgsT2_TRtaGA3 / denomh_sgsT2ta_GA3;
KHTH_STGhT23 = -nomh_sgsT2ta_STGA3 / denomh_sgsT2ta_GA3;
KHTH_GhMOT23 = -nomh_sgsMOT2ta_GA3 / denomh_sgsMOT2ta_GA3;
KHTH_TRGhMOT23 = -nomh_sgsMOT2ta_TRGA3 / denomh_sgsMOT2ta_GA3;
KHTH_TRtaGhMOT23 = -nomh_sgsMOT2_TRtaGA3 / denomh_sgsMOT2ta_GA3;
KHTH_STGhMOT23 = -nomh_sgsMOT2ta_STGA3 / denomh_sgsMOT2ta_GA3;

VEL_sgsta_GA3 = (VEL_sgsta_GA3 / sum_glob3) / AR_fxy_GA;
VEL_sgsta_TRGA3 = (VEL_sgsta_TRGA3 / sum_glob3) / AR_fxy_GA;
VEL_sgsta_STGA3 = (VEL_sgsta_STGA3 / sum_glob3) / AR_fxy_GA;

Lmix_ka3 = KHTH_GhT23 / VEL_sgsta_GA3;
Lmix_kaTR3 = KHTH_TRtaGhT23 / VEL_sgsta_TRGA3;
Lmix_kaST3 = KHTH_STGhT23 / VEL_sgsta_STGA3;
% Lmix3     = -nomh_sgsT2ta_GA3 / denom_MixL_sgsT2ta_GA3;
% Lmix_TR3     = -nomh_sgsT2_TRtaGA3 / denom_MixL_sgsT2_TRtaGA3;
% Lmix_ST3     = -nomh_sgsT2ta_STGA3  / denom_MixL_sgsT2ta_STGA3;

Lmix_MOka3 = KHTH_GhMOT23 / VEL_sgsta_GA3;
Lmix_MOkaTR3 = KHTH_TRtaGhMOT23 / VEL_sgsta_TRGA3;
Lmix_MOkaST3 = KHTH_STGhMOT23 / VEL_sgsta_STGA3;
% Lmix_MO3     = -nomh_sgsMOT2ta_GA3 / denom_MixL_sgsMOT2ta_GA3;
% Lmix_MOTR3     = -nomh_sgsMOT2_TRtaGA3 / denom_MixL_sgsMOT2_TRtaGA3;
% Lmix_MOST3     = -nomh_sgsMOT2ta_STGA3  / denom_MixL_sgsMOT2ta_STGA3;

%%%%%%%%%%%%%%%%%% STANDARD DEVIATION FOR REGRESSION                

%%% NORTH
nomh_sgsT2ta_SD1 = nomh_sgsT2ta_G - ( (nomh_sgsT2ta_GA1 / denomh_sgsT2ta_GA1) * denomh_sgsT2ta_G );
nomh_sgsT2_TRta_SD1 = nomh_sgsT2_TRtaG -  ( (nomh_sgsT2_TRtaGA1 / denomh_sgsT2ta_GA1) * denomh_sgsT2ta_G );          
nomh_sgsT2ta_ST_SD1 = nomh_sgsT2ta_STG - ( (nomh_sgsT2ta_STGA1 / denomh_sgsT2ta_GA1) * denomh_sgsT2ta_G );
nomh_sgsMOT2ta_SD1 = nomh_sgsMOT2ta_G - ( (nomh_sgsMOT2ta_GA1 / denomh_sgsMOT2ta_GA1) * denomh_sgsMOT2ta_G );
nomh_sgsMOT2_TRta_SD1 = nomh_sgsMOT2_TRtaG -  ( (nomh_sgsMOT2_TRtaGA1 / denomh_sgsMOT2ta_GA1) * denomh_sgsMOT2ta_G ); 
nomh_sgsMOT2ta_ST_SD1 = nomh_sgsMOT2ta_STG - ( (nomh_sgsMOT2ta_STGA1 / denomh_sgsMOT2ta_GA1) *  denomh_sgsMOT2ta_G);
denomh_sgsT2ta_SD1 = denomh_sgsT2ta_G - ( (denomh_sgsT2ta_GA1 / sum_glob1 )^0.5 * denomh_sgsT2ta_G.^0.5 );             
denomh_sgsMOT2ta_SD1 = denomh_sgsMOT2ta_G - ( (denomh_sgsMOT2ta_GA1 / sum_glob1)^0.5 * denomh_sgsMOT2ta_G.^0.5 );

nom_MixL_sgsT2ta_SD1 = nom_MixL_sgsT2ta_G - ( (nom_MixL_sgsT2ta_GA1 / denom_MixL_sgsT2ta_GA1) * denom_MixL_sgsT2ta_G );
nom_MixL_sgsT2_TRta_SD1 = nom_MixL_sgsT2_TRtaG -  ( (nom_MixL_sgsT2_TRtaGA1 / denom_MixL_sgsT2ta_GA1) * denom_MixL_sgsT2ta_G );          
nom_MixL_sgsT2ta_ST_SD1 = nom_MixL_sgsT2ta_STG - ( (nom_MixL_sgsT2ta_STGA1 / denom_MixL_sgsT2ta_GA1) * denom_MixL_sgsT2ta_G );
nom_MixL_sgsMOT2ta_SD1 = nom_MixL_sgsMOT2ta_G - ( (nom_MixL_sgsMOT2ta_GA1 / denom_MixL_sgsMOT2ta_GA1) * denom_MixL_sgsMOT2ta_G );
nom_MixL_sgsMOT2_TRta_SD1 = nom_MixL_sgsMOT2_TRtaG -  ( (nom_MixL_sgsMOT2_TRtaGA1 / denom_MixL_sgsMOT2ta_GA1) * denom_MixL_sgsMOT2ta_G ); 
nom_MixL_sgsMOT2ta_ST_SD1 = nom_MixL_sgsMOT2ta_STG - ( (nom_MixL_sgsMOT2ta_STGA1 / denom_MixL_sgsMOT2ta_GA1) *  denom_MixL_sgsMOT2ta_G);
denom_MixL_sgsT2ta_SD1 = denom_MixL_sgsT2ta_G - ( (denom_MixL_sgsT2ta_GA1 / sum_glob1 )^0.5 * denom_MixL_sgsT2ta_G.^0.5 );             
denom_MixL_sgsMOT2ta_SD1 = denom_MixL_sgsMOT2ta_G - ( (denom_MixL_sgsMOT2ta_GA1 / sum_glob1)^0.5 * denom_MixL_sgsMOT2ta_G.^0.5 );

nomh_sgsT2ta_SD1 = nomh_sgsT2ta_SD1 .* nomh_sgsT2ta_SD1 .* MASK_GA .* MASK;
nomh_sgsT2_TRta_SD1 = nomh_sgsT2_TRta_SD1 .* nomh_sgsT2_TRta_SD1 .* MASK_GA .* MASK;
nomh_sgsT2ta_ST_SD1 = nomh_sgsT2ta_ST_SD1 .* nomh_sgsT2ta_ST_SD1 .* MASK_GA .* MASK;
nomh_sgsMOT2ta_SD1 = nomh_sgsMOT2ta_SD1 .* nomh_sgsMOT2ta_SD1 .* MASK_GA .* MASK;
nomh_sgsMOT2_TRta_SD1 = nomh_sgsMOT2_TRta_SD1 .* nomh_sgsMOT2_TRta_SD1 .* MASK_GA .* MASK;
nomh_sgsMOT2ta_ST_SD1 = nomh_sgsMOT2ta_ST_SD1 .* nomh_sgsMOT2ta_ST_SD1 .* MASK_GA .* MASK;
denomh_sgsT2ta_SD1 = denomh_sgsT2ta_SD1 .* denomh_sgsT2ta_SD1 .* MASK_GA .* MASK;
denomh_sgsMOT2ta_SD1 = denomh_sgsMOT2ta_SD1 .* denomh_sgsMOT2ta_SD1 .* MASK_GA .* MASK;

nom_MixL_sgsT2ta_SD1 = nom_MixL_sgsT2ta_SD1 .* nom_MixL_sgsT2ta_SD1 .* MASK_GA .* MASK ;
nom_MixL_sgsT2_TRta_SD1 = nom_MixL_sgsT2_TRta_SD1 .* nom_MixL_sgsT2_TRta_SD1 .* MASK_GA .* MASK;
nom_MixL_sgsT2ta_ST_SD1 = nom_MixL_sgsT2ta_ST_SD1 .* nom_MixL_sgsT2ta_ST_SD1 .* MASK_GA .* MASK;
nom_MixL_sgsMOT2ta_SD1 = nom_MixL_sgsMOT2ta_SD1 .* nom_MixL_sgsMOT2ta_SD1 .* MASK_GA .* MASK;
nom_MixL_sgsMOT2_TRta_SD1 = nom_MixL_sgsMOT2_TRta_SD1 .* nom_MixL_sgsMOT2_TRta_SD1 .* MASK_GA .* MASK;
nom_MixL_sgsMOT2ta_ST_SD1 = nom_MixL_sgsMOT2ta_ST_SD1 .* nom_MixL_sgsMOT2ta_ST_SD1 .* MASK_GA .* MASK;
denom_MixL_sgsT2ta_SD1 = denom_MixL_sgsT2ta_SD1 .* denom_MixL_sgsT2ta_SD1 .* MASK_GA .* MASK;
denom_MixL_sgsMOT2ta_SD1 = denom_MixL_sgsMOT2ta_SD1 .* denom_MixL_sgsMOT2ta_SD1 .* MASK_GA .* MASK;

%%% SOUTH
nomh_sgsT2ta_SD2 = nomh_sgsT2ta_G - ( (nomh_sgsT2ta_GA2 / denomh_sgsT2ta_GA2) * denomh_sgsT2ta_G );
nomh_sgsT2_TRta_SD2 = nomh_sgsT2_TRtaG -  ( (nomh_sgsT2_TRtaGA2 / denomh_sgsT2ta_GA2) * denomh_sgsT2ta_G );          
nomh_sgsT2ta_ST_SD2 = nomh_sgsT2ta_STG - ( (nomh_sgsT2ta_STGA2 / denomh_sgsT2ta_GA2) * denomh_sgsT2ta_G );
nomh_sgsMOT2ta_SD2 = nomh_sgsMOT2ta_G - ( (nomh_sgsMOT2ta_GA2 / denomh_sgsMOT2ta_GA2) * denomh_sgsMOT2ta_G );
nomh_sgsMOT2_TRta_SD2 = nomh_sgsMOT2_TRtaG -  ( (nomh_sgsMOT2_TRtaGA2 / denomh_sgsMOT2ta_GA2) * denomh_sgsMOT2ta_G ); 
nomh_sgsMOT2ta_ST_SD2 = nomh_sgsMOT2ta_STG - ( (nomh_sgsMOT2ta_STGA2 / denomh_sgsMOT2ta_GA2) *  denomh_sgsMOT2ta_G);
denomh_sgsT2ta_SD2 = denomh_sgsT2ta_G - ( (denomh_sgsT2ta_GA2 / sum_glob2 )^0.5 * denomh_sgsT2ta_G.^0.5 );             
denomh_sgsMOT2ta_SD2 = denomh_sgsMOT2ta_G - ( (denomh_sgsMOT2ta_GA2 / sum_glob2)^0.5 * denomh_sgsMOT2ta_G.^0.5 );

nom_MixL_sgsT2ta_SD2 = nom_MixL_sgsT2ta_G - ( (nom_MixL_sgsT2ta_GA2 / denom_MixL_sgsT2ta_GA2) * denom_MixL_sgsT2ta_G );
nom_MixL_sgsT2_TRta_SD2 = nom_MixL_sgsT2_TRtaG -  ( (nom_MixL_sgsT2_TRtaGA2 / denom_MixL_sgsT2ta_GA2) * denom_MixL_sgsT2ta_G );          
nom_MixL_sgsT2ta_ST_SD2 = nom_MixL_sgsT2ta_STG - ( (nom_MixL_sgsT2ta_STGA2 / denom_MixL_sgsT2ta_GA2) * denom_MixL_sgsT2ta_G );
nom_MixL_sgsMOT2ta_SD2 = nom_MixL_sgsMOT2ta_G - ( (nom_MixL_sgsMOT2ta_GA2 / denom_MixL_sgsMOT2ta_GA2) * denom_MixL_sgsMOT2ta_G );
nom_MixL_sgsMOT2_TRta_SD2 = nom_MixL_sgsMOT2_TRtaG -  ( (nom_MixL_sgsMOT2_TRtaGA2 / denom_MixL_sgsMOT2ta_GA2) * denom_MixL_sgsMOT2ta_G ); 
nom_MixL_sgsMOT2ta_ST_SD2 = nom_MixL_sgsMOT2ta_STG - ( (nom_MixL_sgsMOT2ta_STGA2 / denom_MixL_sgsMOT2ta_GA2) *  denom_MixL_sgsMOT2ta_G);
denom_MixL_sgsT2ta_SD2 = denom_MixL_sgsT2ta_G - ( (denom_MixL_sgsT2ta_GA2 / sum_glob2 )^0.5 * denom_MixL_sgsT2ta_G.^0.5 );             
denom_MixL_sgsMOT2ta_SD2 = denom_MixL_sgsMOT2ta_G - ( (denom_MixL_sgsMOT2ta_GA2 / sum_glob2)^0.5 * denom_MixL_sgsMOT2ta_G.^0.5 );

nomh_sgsT2ta_SD2 = nomh_sgsT2ta_SD2 .* nomh_sgsT2ta_SD2 .* MASK_GA .* MASK ;
nomh_sgsT2_TRta_SD2 = nomh_sgsT2_TRta_SD2 .* nomh_sgsT2_TRta_SD2 .* MASK_GA .* MASK;
nomh_sgsT2ta_ST_SD2 = nomh_sgsT2ta_ST_SD2 .* nomh_sgsT2ta_ST_SD2 .* MASK_GA .* MASK;
nomh_sgsMOT2ta_SD2 = nomh_sgsMOT2ta_SD2 .* nomh_sgsMOT2ta_SD2 .* MASK_GA .* MASK;
nomh_sgsMOT2_TRta_SD2 = nomh_sgsMOT2_TRta_SD2 .* nomh_sgsMOT2_TRta_SD2 .* MASK_GA .* MASK;
nomh_sgsMOT2ta_ST_SD2 = nomh_sgsMOT2ta_ST_SD2 .* nomh_sgsMOT2ta_ST_SD2 .* MASK_GA .* MASK;
denomh_sgsT2ta_SD2 = denomh_sgsT2ta_SD2 .* denomh_sgsT2ta_SD2 .* MASK_GA .* MASK;
denomh_sgsMOT2ta_SD2 = denomh_sgsMOT2ta_SD2 .* denomh_sgsMOT2ta_SD2 .* MASK_GA .* MASK;

nom_MixL_sgsT2ta_SD2 = nom_MixL_sgsT2ta_SD2 .* nom_MixL_sgsT2ta_SD2 .* MASK_GA .* MASK ;
nom_MixL_sgsT2_TRta_SD2 = nom_MixL_sgsT2_TRta_SD2 .* nom_MixL_sgsT2_TRta_SD2 .* MASK_GA .* MASK;
nom_MixL_sgsT2ta_ST_SD2 = nom_MixL_sgsT2ta_ST_SD2 .* nom_MixL_sgsT2ta_ST_SD2 .* MASK_GA .* MASK;
nom_MixL_sgsMOT2ta_SD2 = nom_MixL_sgsMOT2ta_SD2 .* nom_MixL_sgsMOT2ta_SD2 .* MASK_GA .* MASK;
nom_MixL_sgsMOT2_TRta_SD2 = nom_MixL_sgsMOT2_TRta_SD2 .* nom_MixL_sgsMOT2_TRta_SD2 .* MASK_GA .* MASK;
nom_MixL_sgsMOT2ta_ST_SD2 = nom_MixL_sgsMOT2ta_ST_SD2 .* nom_MixL_sgsMOT2ta_ST_SD2 .* MASK_GA .* MASK;
denom_MixL_sgsT2ta_SD2 = denom_MixL_sgsT2ta_SD2 .* denom_MixL_sgsT2ta_SD2 .* MASK_GA .* MASK;
denom_MixL_sgsMOT2ta_SD2 = denom_MixL_sgsMOT2ta_SD2 .* denom_MixL_sgsMOT2ta_SD2 .* MASK_GA .* MASK;



nomh_sgsT2ta_SD1_GA = 0;
nomh_sgsT2_TRta_SD1_GA = 0;
nomh_sgsT2ta_ST_SD1_GA = 0;
nomh_sgsMOT2ta_SD1_GA = 0;
nomh_sgsMOT2_TRta_SD1_GA = 0;              
nomh_sgsMOT2ta_ST_SD1_GA  = 0;
denomh_sgsT2ta_SD1_GA = 0;      
denomh_sgsMOT2ta_SD1_GA = 0;

nom_MixL_sgsT2ta_SD1_GA = 0;
nom_MixL_sgsT2_TRta_SD1_GA = 0;
nom_MixL_sgsT2ta_ST_SD1_GA = 0;
nom_MixL_sgsMOT2ta_SD1_GA = 0;
nom_MixL_sgsMOT2_TRta_SD1_GA = 0;              
nom_MixL_sgsMOT2ta_ST_SD1_GA  = 0;
denom_MixL_sgsT2ta_SD1_GA = 0;      
denom_MixL_sgsMOT2ta_SD1_GA = 0;

nomh_sgsT2ta_SD2_GA = 0;
nomh_sgsT2_TRta_SD2_GA = 0;
nomh_sgsT2ta_ST_SD2_GA = 0;
nomh_sgsMOT2ta_SD2_GA = 0;
nomh_sgsMOT2_TRta_SD2_GA = 0;              
nomh_sgsMOT2ta_ST_SD2_GA  = 0;
denomh_sgsT2ta_SD2_GA = 0;      
denomh_sgsMOT2ta_SD2_GA = 0;

nom_MixL_sgsT2ta_SD2_GA = 0;
nom_MixL_sgsT2_TRta_SD2_GA = 0;
nom_MixL_sgsT2ta_ST_SD2_GA = 0;
nom_MixL_sgsMOT2ta_SD2_GA = 0;
nom_MixL_sgsMOT2_TRta_SD2_GA = 0;              
nom_MixL_sgsMOT2ta_ST_SD2_GA  = 0;
denom_MixL_sgsT2ta_SD2_GA = 0;      
denom_MixL_sgsMOT2ta_SD2_GA = 0;
for k=1:zsize
    for i=1:(xsize_f)  %%%
        for j=1:(ysize_f)
            
             if(Y_f(j) >= -33.3)    %%% -33.3 For Neverland
                 nomh_sgsT2ta_SD1_GA = nomh_sgsT2ta_SD1_GA + nomh_sgsT2ta_SD1(i,j,k);
                 nomh_sgsT2_TRta_SD1_GA = nomh_sgsT2_TRta_SD1_GA + nomh_sgsT2_TRta_SD1(i,j,k);
                 nomh_sgsT2ta_ST_SD1_GA = nomh_sgsT2ta_ST_SD1_GA + nomh_sgsT2ta_ST_SD1(i,j,k);
                 nomh_sgsMOT2ta_SD1_GA =  nomh_sgsMOT2ta_SD1_GA  + nomh_sgsMOT2ta_SD1(i,j,k);
                 nomh_sgsMOT2_TRta_SD1_GA = nomh_sgsMOT2_TRta_SD1_GA + nomh_sgsMOT2_TRta_SD1(i,j,k);
                 nomh_sgsMOT2ta_ST_SD1_GA = nomh_sgsMOT2ta_ST_SD1_GA + nomh_sgsMOT2ta_ST_SD1(i,j,k);
                 denomh_sgsT2ta_SD1_GA = denomh_sgsT2ta_SD1_GA + denomh_sgsT2ta_SD1(i,j,k);
                 denomh_sgsMOT2ta_SD1_GA = denomh_sgsMOT2ta_SD1_GA  + denomh_sgsMOT2ta_SD1(i,j,k);
                 
                  
                 nom_MixL_sgsT2ta_SD1_GA = nom_MixL_sgsT2ta_SD1_GA + nom_MixL_sgsT2ta_SD1(i,j,k);
                 nom_MixL_sgsT2_TRta_SD1_GA = nom_MixL_sgsT2_TRta_SD1_GA + nom_MixL_sgsT2_TRta_SD1(i,j,k);
                 nom_MixL_sgsT2ta_ST_SD1_GA = nom_MixL_sgsT2ta_ST_SD1_GA + nom_MixL_sgsT2ta_ST_SD1(i,j,k);
                 nom_MixL_sgsMOT2ta_SD1_GA =  nom_MixL_sgsMOT2ta_SD1_GA  + nom_MixL_sgsMOT2ta_SD1(i,j,k);
                 nom_MixL_sgsMOT2_TRta_SD1_GA = nom_MixL_sgsMOT2_TRta_SD1_GA + nom_MixL_sgsMOT2_TRta_SD1(i,j,k);
                 nom_MixL_sgsMOT2ta_ST_SD1_GA = nom_MixL_sgsMOT2ta_ST_SD1_GA + nom_MixL_sgsMOT2ta_ST_SD1(i,j,k);
                 denom_MixL_sgsT2ta_SD1_GA = denom_MixL_sgsT2ta_SD1_GA + denom_MixL_sgsT2ta_SD1(i,j,k);
                 denom_MixL_sgsMOT2ta_SD1_GA = denom_MixL_sgsMOT2ta_SD1_GA  + denom_MixL_sgsMOT2ta_SD1(i,j,k);
            
             elseif (Y_f(j) < -50)  %%% -50 For Neverland 
                 
                 nomh_sgsT2ta_SD2_GA = nomh_sgsT2ta_SD2_GA + nomh_sgsT2ta_SD2(i,j,k);
                 nomh_sgsT2_TRta_SD2_GA = nomh_sgsT2_TRta_SD2_GA + nomh_sgsT2_TRta_SD2(i,j,k);
                 nomh_sgsT2ta_ST_SD2_GA = nomh_sgsT2ta_ST_SD2_GA + nomh_sgsT2ta_ST_SD2(i,j,k);
                 nomh_sgsMOT2ta_SD2_GA =  nomh_sgsMOT2ta_SD2_GA  + nomh_sgsMOT2ta_SD2(i,j,k);
                 nomh_sgsMOT2_TRta_SD2_GA = nomh_sgsMOT2_TRta_SD2_GA + nomh_sgsMOT2_TRta_SD2(i,j,k);
                 nomh_sgsMOT2ta_ST_SD2_GA = nomh_sgsMOT2ta_ST_SD2_GA + nomh_sgsMOT2ta_ST_SD2(i,j,k);
                 denomh_sgsT2ta_SD2_GA = denomh_sgsT2ta_SD2_GA + denomh_sgsT2ta_SD2(i,j,k);
                 denomh_sgsMOT2ta_SD2_GA = denomh_sgsMOT2ta_SD2_GA  + denomh_sgsMOT2ta_SD2(i,j,k);
                             
                 
                 nom_MixL_sgsT2ta_SD2_GA = nom_MixL_sgsT2ta_SD2_GA + nom_MixL_sgsT2ta_SD2(i,j,k);
                 nom_MixL_sgsT2_TRta_SD2_GA = nom_MixL_sgsT2_TRta_SD2_GA + nom_MixL_sgsT2_TRta_SD2(i,j,k);
                 nom_MixL_sgsT2ta_ST_SD2_GA = nom_MixL_sgsT2ta_ST_SD2_GA + nom_MixL_sgsT2ta_ST_SD2(i,j,k);
                 nom_MixL_sgsMOT2ta_SD2_GA =  nom_MixL_sgsMOT2ta_SD2_GA  + nom_MixL_sgsMOT2ta_SD2(i,j,k);
                 nom_MixL_sgsMOT2_TRta_SD2_GA = nom_MixL_sgsMOT2_TRta_SD2_GA + nom_MixL_sgsMOT2_TRta_SD2(i,j,k);
                 nom_MixL_sgsMOT2ta_ST_SD2_GA = nom_MixL_sgsMOT2ta_ST_SD2_GA + nom_MixL_sgsMOT2ta_ST_SD2(i,j,k);
                 denom_MixL_sgsT2ta_SD2_GA = denom_MixL_sgsT2ta_SD2_GA + denom_MixL_sgsT2ta_SD2(i,j,k);
                 denom_MixL_sgsMOT2ta_SD2_GA = denom_MixL_sgsMOT2ta_SD2_GA  + denom_MixL_sgsMOT2ta_SD2(i,j,k);
             end
        end
    end
end


KHTH_GhT21_SD =   ( (nomh_sgsT2ta_SD1_GA/sum_glob1) /   denomh_sgsT2ta_SD1_GA )^0.5;     
KHTH_TRtaGhT21_SD = ( (nomh_sgsT2_TRta_SD1_GA/sum_glob1) / denomh_sgsT2ta_SD1_GA )^0.5;   
KHTH_STGhT21_SD = ( (nomh_sgsT2ta_ST_SD1_GA/sum_glob1) / denomh_sgsT2ta_SD1_GA )^0.5;     
KHTH_GhMOT21_SD = ( (nomh_sgsMOT2ta_SD1_GA/sum_glob1) /   denomh_sgsMOT2ta_SD1_GA )^0.5; 
KHTH_TRtaGhMOT21_SD = ( (nomh_sgsMOT2_TRta_SD1_GA/sum_glob1) / denomh_sgsMOT2ta_SD1_GA )^0.5;
KHTH_STGhMOT21_SD = ( (nomh_sgsMOT2ta_ST_SD1_GA/sum_glob1) / denomh_sgsMOT2ta_SD1_GA )^0.5; 

Lmix_T21_SD =   ( (nom_MixL_sgsT2ta_SD1_GA/sum_glob1) /   denom_MixL_sgsT2ta_SD1_GA )^0.5;     
Lmix_TRtaT21_SD = ( (nom_MixL_sgsT2_TRta_SD1_GA/sum_glob1) / denom_MixL_sgsT2ta_SD1_GA )^0.5;   
Lmix_STT21_SD = ( (nom_MixL_sgsT2ta_ST_SD1_GA/sum_glob1) / denom_MixL_sgsT2ta_SD1_GA )^0.5;     
Lmix_MOT21_SD = ( (nom_MixL_sgsMOT2ta_SD1_GA/sum_glob1) /   denom_MixL_sgsMOT2ta_SD1_GA )^0.5; 
Lmix_TRtaMOT21_SD = ( (nom_MixL_sgsMOT2_TRta_SD1_GA/sum_glob1) / denom_MixL_sgsMOT2ta_SD1_GA )^0.5;
Lmix_STMOT21_SD = ( (nom_MixL_sgsMOT2ta_ST_SD1_GA/sum_glob1) / denom_MixL_sgsMOT2ta_SD1_GA )^0.5; 

KHTH_GhT22_SD =   ( (nomh_sgsT2ta_SD2_GA/sum_glob2) /   denomh_sgsT2ta_SD2_GA )^0.5;     
KHTH_TRtaGhT22_SD = ( (nomh_sgsT2_TRta_SD2_GA/sum_glob2) / denomh_sgsT2ta_SD2_GA )^0.5;   
KHTH_STGhT22_SD = ( (nomh_sgsT2ta_ST_SD2_GA/sum_glob2) / denomh_sgsT2ta_SD2_GA )^0.5;     
KHTH_GhMOT22_SD = ( (nomh_sgsMOT2ta_SD2_GA/sum_glob2) /   denomh_sgsMOT2ta_SD2_GA )^0.5; 
KHTH_TRtaGhMOT22_SD = ( (nomh_sgsMOT2_TRta_SD2_GA/sum_glob2) / denomh_sgsMOT2ta_SD2_GA )^0.5;
KHTH_STGhMOT22_SD = ( (nomh_sgsMOT2ta_ST_SD2_GA/sum_glob2) / denomh_sgsMOT2ta_SD2_GA )^0.5;

Lmix_T22_SD =   ( (nom_MixL_sgsT2ta_SD2_GA/sum_glob2) /   denom_MixL_sgsT2ta_SD2_GA )^0.5;     
Lmix_TRtaT22_SD = ( (nom_MixL_sgsT2_TRta_SD2_GA/sum_glob2) / denom_MixL_sgsT2ta_SD2_GA )^0.5;   
Lmix_STT22_SD = ( (nom_MixL_sgsT2ta_ST_SD2_GA/sum_glob2) / denom_MixL_sgsT2ta_SD2_GA )^0.5;     
Lmix_MOT22_SD = ( (nom_MixL_sgsMOT2ta_SD2_GA/sum_glob2) /   denom_MixL_sgsMOT2ta_SD2_GA )^0.5; 
Lmix_TRtaMOT22_SD = ( (nom_MixL_sgsMOT2_TRta_SD2_GA/sum_glob2) / denom_MixL_sgsMOT2ta_SD2_GA )^0.5;
Lmix_STMOT22_SD = ( (nom_MixL_sgsMOT2ta_ST_SD2_GA/sum_glob2) / denom_MixL_sgsMOT2ta_SD2_GA )^0.5; 

%%%%%% 95% confidence 
z_s=1.96;
KHTH_GhT21_p95 = z_s * KHTH_GhT21_SD;
KHTH_TRtaGhT21_p95 = z_s * KHTH_TRtaGhT21_SD;
KHTH_STGhT21_p95 = z_s * KHTH_STGhT21_SD;
KHTH_GhMOT21_p95 = z_s * KHTH_GhMOT21_SD;
KHTH_TRtaGhMOT21_p95 = z_s * KHTH_TRtaGhMOT21_SD;
KHTH_STGhMOT21_p95 = z_s * KHTH_STGhMOT21_SD;

Lmix_T21_p95 = z_s * Lmix_T21_SD;
Lmix_TRtaT21_p95 = z_s * Lmix_TRtaT21_SD;
Lmix_STT21_p95 = z_s * Lmix_STT21_SD;
Lmix_MOT21_p95 = z_s * Lmix_MOT21_SD;
Lmix_TRtaMOT21_p95 = z_s * Lmix_TRtaMOT21_SD;
Lmix_STMOT21_p95 = z_s * Lmix_STMOT21_SD;

KHTH_GhT22_p95 = z_s * KHTH_GhT22_SD;
KHTH_TRtaGhT22_p95 = z_s * KHTH_TRtaGhT22_SD;
KHTH_STGhT22_p95 = z_s * KHTH_STGhT22_SD;
KHTH_GhMOT22_p95 = z_s * KHTH_GhMOT22_SD;
KHTH_TRtaGhMOT22_p95 = z_s * KHTH_TRtaGhMOT22_SD;
KHTH_STGhMOT22_p95 = z_s * KHTH_STGhMOT22_SD;

Lmix_T22_p95 = z_s * Lmix_T22_SD;
Lmix_TRtaT22_p95 = z_s * Lmix_TRtaT22_SD;
Lmix_STT22_p95 = z_s * Lmix_STT22_SD;
Lmix_MOT22_p95 = z_s * Lmix_MOT22_SD;
Lmix_TRtaMOT22_p95 = z_s * Lmix_TRtaMOT22_SD;
Lmix_STMOT22_p95 = z_s * Lmix_STMOT22_SD;





toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MASK_CORP_NWR = ones(ysize_f,zsize,'single');  %%%%%% MAKS CORUPTION for NWR
for j=1:ysize_f
    if (Y_f(j)<-52.00)
        MASK_CORP_NWR (j,1) = 0;
    end
    if (Y_f(j)<-58.00)
        MASK_CORP_NWR (j,2) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return


%%%%%%%%%% Deformation radius
figure(4345)
contourf(X,Y',RD1_ta'/1000,[5:1:55]);
set(gca,'clim',[5,55]);
colorbar;
hold on

figure(23086) %%% Deformation radius
plot(Y,RD1_taza/1000,'-r','linewidth',2);
hold on



