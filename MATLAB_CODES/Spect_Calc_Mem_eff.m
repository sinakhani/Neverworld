%%%%%%%% THIS SCRIPT GENERATES KINETIC (POTENTIAL) ENERGY, EDDY KINETIC ENERGY and SUBGRID SCALE EKE IN THE CHANNEL AND NEVERWORLD 
%%%%%%%% CONFIGRATIONS USING MOM6 NETCDF DATA. PLEASE SEE THE PAPER:
%%%%%%%% (https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019MS001721)
%%%%%%%% THIS SCRIPT IS WRITTEN BY SINA KHANI (skhani@uw.edu).
%%%%%%%% THIS CODE IS FREE TO USE SUBJECT TO PROPER CITATION OF THE WORK.

openstr = 'ocean_geometry.nc'; %%%%% READ OCEAN GEOMETRY NETCDF FILE HERE.
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
Mean_DY = mean(DY,1)';
DYY = DY(1,:)';

openstr = 'snapshots.nc'; %%%% READ SNAPSHOT NETCDF FILE HERE.
fi = netcdf.open(openstr, 'nowrite');
% get the varibles to plot, there are many to choose from
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

AR(AR==0) = nan;
AR_GA =0;
sum_areaav = 0;
for i=1:xsize
    for j=1:ysize
        if (~isnan(AR(i,j)))
            sum_areaav = sum_areaav + 1;
            AR_GA = AR_GA + AR(i,j);
        end
    end
end
AR_GA = AR_GA / sum_areaav;
AR(isnan(AR)) = 0;
ARR = AR(1,:)';


tsize = 730;%%1825; %%730 Maximum number of years for <= 1/4 deg tsize=3650 and for 1/16 tsize=730, tsize=50 for 1/8 deg. 730 for 1/16 %%% 1825 for 1/8
nyr = 10;   %% Do not change this at all! 
avyr = 0.25;  %% !!!!!! for 1/8 in movmean use 0.25 and for 1/16 and 1/4 (3650) use 0.25/4 
TS = 365;%%1460;%%365;%floor( tsize * ((nyr-avyr)/nyr) ); %% time average for TS<t<tsize

gprime = [9.81,0.015,0.005,0.003,0.002,0.001];


X=netcdf.getVar(fi,0);
Y=netcdf.getVar(fi,1);
Z=netcdf.getVar(fi,2);
%lon = zeros(xsize,ysize,'single');


y = Y * (111.32 * 1000); %[m]
x = X * (111.32 * 1000); % [m]
z = [0 -150 -400 -1000 -2000 -3000];
%z = fliplr(z);
thcknss = [150 250 600 1000 1000 1000];
H=sum(thcknss);
thcknss = thcknss/H;
OMEGA = 7.2921e-5;

% 
% for i=1:xsize
%     for j=1:ysize
%         lon(i,j) = x(i) * cosd(Y(j));
%     end
% end

  
    
lon = x(xsize) * cosd(50);
lon1 = x(xsize) * cosd(23.75);
lon2 = x(xsize) * cosd(49.15);


CORCOEF = zeros(xsize,ysize,1,'single');
CORCOEF1 = zeros(ysize,1,'single');
beta = zeros(ysize,1,'single');


for j=1:ysize
    CORCOEF(:,j) = 2 * OMEGA * sind(Y(j));  %2 * OMEGA * sind(Y(j)); % -4.985E-5 for Y=-20 % -1.1172E-4 for Y=-50;
    CORCOEF1(j) = 2 * OMEGA * sind(Y(j));
end



u = zeros(xsize,ysize,zsize,1,'single');
v = zeros(xsize,ysize,zsize,1,'single');
h = zeros(xsize,ysize,zsize,1,'single');
e = zeros(xsize,ysize,zsize,1,'single');

uta = zeros(xsize,ysize,zsize,'single');
vta = zeros(xsize,ysize,zsize,'single');
ugta = zeros(xsize,ysize,zsize,'single');
vgta = zeros(xsize,ysize,zsize,'single');
hta = zeros(xsize,ysize,zsize,'single');
eta = zeros(xsize,ysize,zsize,'single');
Enr_realta = zeros(xsize,ysize,zsize,'single');
VEL_realta = zeros(xsize,ysize,zsize,'single');
Enr_realta_ED = zeros(xsize,ysize,zsize,'single');
VEL_realta_ED = zeros(xsize,ysize,zsize,'single');
VEL_realta_ED_BAR = zeros(xsize,ysize,zsize,'single');
VEL_realta_BAR = zeros(xsize,ysize,zsize,'single');
U_tta = zeros(xsize,ysize,'single');
V_tta = zeros(xsize,ysize,'single');
VEL_Tta =  zeros(xsize,ysize,'single');
VEL_Bta =  zeros(xsize,ysize,'single');


hta_za= zeros(xsize,ysize,'single');

uprime = zeros(xsize,ysize,zsize,'single');
vprime = zeros(xsize,ysize,zsize,'single');
hprime = zeros(xsize,ysize,zsize,'single');
eprime = zeros(xsize,ysize,zsize,'single');
ugprime = zeros(xsize,ysize,zsize,'single');
vgprime = zeros(xsize,ysize,zsize,'single');

Enr_EDta  = zeros(xsize,ysize,zsize,'single');
PEnr_EDta = zeros(xsize,ysize,zsize,'single');

Enrta  = zeros(xsize,ysize,zsize,'single');
PEnrta = zeros(xsize,ysize,zsize,'single');

ENG_ed = zeros(xsize,ysize,zsize,'single');

for m=1:ysize
                if (m==1)
                    beta(m) = (CORCOEF1(m+1) - CORCOEF1(m)) / (DY(m));
                elseif (m>1) && (m<ysize)
                    beta(m) = (CORCOEF1(m+1) - CORCOEF1(m-1)) / (2*DY(m));
                elseif (m==ysize)
                    beta(m) = (CORCOEF1(m) - CORCOEF1(m-1)) / (DY(m)); 
                end                
end


tic
sum_a = 0;

MASK = ones(size(h));
for i =1:xsize
    for j=1:ysize
        for k=1:zsize
            if (k==1 && Y(j) < -40.22)
                MASK(i,j,k) =0;
            elseif (k==2 && Y(j) < -44.28)
                MASK(i,j,k) =0;
            elseif (k==3 && Y(j) < -52.48)
                MASK(i,j,k) =0;
            end
        end
    end
end

for i= 1:tsize    
    if (i>=TS)   
        u = ncread(openstr,'u',[1,1,1,i],[xsize,ysize,zsize,1]);    
        v = ncread(openstr,'v',[1,1,1,i],[xsize,ysize,zsize,1]);    
        h = ncread(openstr,'h',[1,1,1,i],[xsize,ysize,zsize,1]);     
        e = ncread(openstr,'e',[1,1,1,i],[xsize,ysize,zsize,1]);
        eee2 = ncread(openstr,'e',[1,1,1,i],[xsize,ysize,zsize+1,1]);
        
        u(u<-1.e33) = NaN;
        v(v<-1.e33) = NaN;
        h(h<-1.e33) = NaN;
        h(h==0.0) = NaN;
        e(e<-1.e33) = NaN;
        eee2(eee2<-1.e33) = NaN;
        [ug,vg]=spc_geovel(eee2,CORCOEF,gprime,DX,DY,eps,1);

        
        
            
        sum_a = sum_a +1;
        uta(:,:,:) = uta(:,:,:) + u(:,:,:);
        vta(:,:,:) = vta(:,:,:) + v(:,:,:);
        hta(:,:,:) = hta(:,:,:) + h(:,:,:);
        eta(:,:,:) = eta(:,:,:) + e(:,:,:);
        ugta(:,:,:) = ugta(:,:,:) + ug(:,:,:);
        vgta(:,:,:) = vgta(:,:,:) + vg(:,:,:);
    end
end


uta = uta / (sum_a);
vta = vta / (sum_a);
hta = hta / (sum_a);
eta = eta / (sum_a);
ugta = ugta / (sum_a);
vgta = vgta / (sum_a);

mean_hta = nanmean(hta(:));

ugta_za = nanmean(ugta,1);
vgta_za = nanmean(vgta,1);
ugta_za = permute(ugta_za,[2,3,1]);
vgta_za = permute(vgta_za,[2,3,1]);

uta(isnan(uta))=0;
vta(isnan(vta))=0;
hta(isnan(hta))=0;
eta(isnan(eta))=0;
ugta(isnan(ugta))=0;
vgta(isnan(vgta))=0;
ugta_za(isnan(ugta_za))=0;
vgta_za(isnan(vgta_za))=0;




toc

tic
eps=5;
sum_ap = 0;

for i= 1:tsize   
    if (i>=TS)    
        u = ncread(openstr,'u',[1,1,1,i],[xsize,ysize,zsize,1]);    
        v = ncread(openstr,'v',[1,1,1,i],[xsize,ysize,zsize,1]); 
        h = ncread(openstr,'h',[1,1,1,i],[xsize,ysize,zsize,1]); 
        e = ncread(openstr,'e',[1,1,1,i],[xsize,ysize,zsize,1]);  
        eee = ncread(openstr,'e',[1,1,1,i],[xsize,ysize,zsize+1,1]);
            
        uprime =  u - uta;
        vprime =  v - vta;
        eprime =  e - eta;
        
        
        uprime(uprime<-1.e33) = NaN;
        vprime(vprime<-1.e33) = NaN;
        eprime(eprime<-1.e33) = NaN;
        
        u(u<-1.e33) = NaN;
        v(v<-1.e33) = NaN;
        e(e<-1.e33) = NaN;
        
        [ug,vg]=spc_geovel(eee,CORCOEF,gprime,DX,DY,eps,1);
        
        for m=1:ysize
            for n=1:zsize
                ugprime(:,m,n) = ug(:,m,n) - ugta_za(m,n);
                vgprime(:,m,n) = vg(:,m,n) - vgta_za(m,n);
            end
        end
        
        uprime(isnan(uprime))=0;
        vprime(isnan(vprime))=0;
        eprime(isnan(eprime))=0;
        ugprime(isnan(ugprime))=0;
        vgprime(isnan(vgprime))=0;
        
                
        
        u(isnan(u))=0;
        v(isnan(v))=0;
        e(isnan(e))=0;
        ug(isnan(ug))=0;
        vg(isnan(vg))=0;
        
                

        uprimek = fft(ugprime);
        vprimek = fft(vgprime);
        eprimek = fft(eprime);
        
        uk = fft(ug); %fft(ugprime);  %%% Change to geostrophic velocity 
        vk = fft(vg); %fft(vgprime);  %%% Change to geostrophic velocity 
        ek = fft(e);
    
        uprimek = uprimek / xsize;
        vprimek = vprimek / xsize;
        eprimek = eprimek / xsize;
        
        uk = uk / xsize;
        vk = vk / xsize;
        ek = ek / xsize;
        
        Enr_real =  (ug .* ug + vg .* vg );  %%% (ugprime .* ugprime + vgprime .* vgprime );
        VEL_real = Enr_real;
        
        VEL_real_h = (VEL_real .* h) .* MASK;
        ut_h = (ug .* h) .* MASK;
        vt_h = (vg .* h) .* MASK;
        
        sum_VEL_real_h = sum(VEL_real_h,3);
        sum_ut_h = sum(ut_h,3);
        sum_vt_h = sum(vt_h,3);
        sum_h = sum(h,3);
        
        U_t = sum_ut_h ./ sum_h;
        V_t = sum_vt_h ./ sum_h;
        VEL_T = sum_VEL_real_h ./ sum_h; 
        VEL_B = U_t .* U_t + V_t .* V_t;
        
        
        
        
        Enr_real_ED =  (ugprime .* ugprime + vgprime .* vgprime ); %(u .* u + v .* v ); %%% Change to geostrophic velocity
        VEL_real_ED = Enr_real_ED;
        
        
        
    
        Enr_ED   =  ( uprimek  .* conj(uprimek ) + vprimek  .* conj(vprimek ) );
        PEnr_ED =  ( eprimek .* conj(eprimek) );
        
        Enr   = ( uk  .* conj(uk ) + vk  .* conj(vk ) );
        PEnr  =  ( ek .* conj(ek) );
        
        ENG  = ( (uprime .* uprime) + (vprime  .* vprime ) );
        
  
            
        sum_ap = sum_ap +1;
        Enr_EDta(:,:,:)  = Enr_EDta(:,:,:)  + Enr_ED(:,:,:) ;
        PEnr_EDta(:,:,:)  = PEnr_EDta(:,:,:)  + PEnr_ED(:,:,:) ; 
        Enrta(:,:,:)  = Enrta(:,:,:)  + Enr(:,:,:) ;
        PEnrta(:,:,:)  = PEnrta(:,:,:)  + PEnr(:,:,:) ;
     
        
        ENG_ed(:,:,:)  = ENG_ed(:,:,:) + ENG(:,:,:); 
        Enr_realta(:,:,:) =  Enr_realta(:,:,:)  +  Enr_real(:,:,:);
        VEL_realta(:,:,:)  =  VEL_realta(:,:,:)  +  VEL_real(:,:,:);
        Enr_realta_ED(:,:,:) =  Enr_realta_ED(:,:,:)  +  Enr_real_ED(:,:,:);
        VEL_realta_ED(:,:,:) =  VEL_realta_ED(:,:,:)  +  VEL_real_ED(:,:,:);
        U_tta(:,:) =  U_tta(:,:)  +  U_t(:,:);
        V_tta(:,:) =  V_tta(:,:)  +  V_t(:,:);
        VEL_Tta(:,:) =  VEL_Tta(:,:)  +  VEL_T(:,:);
        VEL_Bta(:,:) =  VEL_Bta(:,:)  +  VEL_B(:,:);
        

    end
    
    
end



Enr_EDta  = Enr_EDta / (sum_ap);
PEnr_EDta = PEnr_EDta / (sum_ap);
Enrta  = Enrta / (sum_ap);
PEnrta  = PEnrta / (sum_ap);
ENG_edta = ENG_ed / (sum_ap);

VEL_edta = ENG_edta.^0.5;
lon_m = mean(lon,1)';

Enr_realta = Enr_realta/ (sum_ap);
VEL_realta = VEL_realta/ (sum_ap);

Enr_realta_ED = Enr_realta_ED/ (sum_ap);
VEL_realta_ED = VEL_realta_ED/ (sum_ap);

U_tta = U_tta/ (sum_ap);
V_tta = V_tta/ (sum_ap);
VEL_Tta = VEL_Tta/ (sum_ap);
VEL_Bta = VEL_Bta/ (sum_ap); %U_tta .* U_tta + V_tta .* V_tta;
                


toc


%%%%%%%%%%% Calculating kinetic and potential energy spectrum

utak = fft(uta);
vtak = fft(vta);
etak = fft(eta);

clear uta; clear vta; clear eta;

utak = utak / xsize;
vtak = vtak / xsize;
etak = etak / xsize;

MEnr  = (1/2) * ( utak .* conj(utak) + vtak .* conj(vtak) );
MPEnr = (1/2) * ( etak .* conj(etak) );


bins = xsize;
num = zeros(1,bins);
Enr_sp = zeros(bins,ysize,zsize,'single');
Enr_spint = zeros(bins,ysize,zsize,'single');
PEnr_sp = zeros(bins,ysize,zsize,'single');
EnrED_sp = zeros(bins,ysize,zsize,'single');
PEnrED_sp = zeros(bins,ysize,zsize,'single');
MEnr_sp = zeros(bins,ysize,zsize,'single');
MPEnr_sp = zeros(bins,ysize,zsize,'single');

 
KX = [0:((xsize+1)/2)-1,-((xsize+1)/2):-1];  
kc = floor((xsize+1)/2);

tic
Enrta(Enrta==0) = nan;
for k =1:zsize
    for j=1:ysize
        for i=1:xsize
            if (~isnan(Enrta(i,j,k)))
                wkx = floor(abs(KX(i)));
                p=wkx+1;
                num(p) = p;
                Enr_sp(p,j,k) = Enr_sp(p,j,k) + ( Enrta(i,j,k) * (lon/ (2*pi))*(ARR(j)) );
                Enr_spint(p,j,k) = Enr_spint(p,j,k) + ( Enrta(i,j,k)/ (2*pi) *(ARR(j)) );
                PEnr_sp(p,j,k) = PEnr_sp(p,j,k) + ( PEnrta(i,j,k) * (lon/ (2*pi))*(ARR(j))  );
                EnrED_sp(p,j,k) = EnrED_sp(p,j,k) + ( Enr_EDta(i,j,k)* (lon/ (2*pi))*(ARR(j)) );
                PEnrED_sp(p,j,k) = PEnrED_sp(p,j,k) + ( PEnr_EDta(i,j,k) * (lon / (2*pi)) *(ARR(j)));
                MEnr_sp(p,j,k) = MEnr_sp(p,j,k) + ( MEnr(i,j,k) * (lon/ (2*pi)) *(ARR(j)) );
                MPEnr_sp(p,j,k) = MPEnr_sp(p,j,k) + ( MPEnr(i,j,k) * (lon/ (2*pi)) *(ARR(j)));
            end
        end
    end
end

Enr_sp = Enr_sp / AR_GA;
Enr_spint = Enr_spint / AR_GA;
PEnr_sp = PEnr_sp / AR_GA;
EnrED_sp = EnrED_sp / AR_GA;
PEnrED_sp = PEnrED_sp / AR_GA;
MEnr_sp = MEnr_sp / AR_GA;
MPEnr_sp = MPEnr_sp / AR_GA;


toc

psize = size(num,2);

Enrta(isnan(Enrta)) = 0;


KX_INV = 1./KX;
num_edd_length = zeros(ysize,1);
dnum_edd_length = zeros(ysize,1);
num2_edd_length = zeros(ysize,1);
%%%%%%%%%%%%%%%%% Transient eddy length scale 
for k =1:zsize
    for j=1:ysize
        for i=1:psize
            num_edd_length(j)  = num_edd_length(j)  +  (EnrED_sp(i,j,k) * (lon/ (2*pi)) );
            dnum_edd_length(j) = dnum_edd_length(j) +  ( KX(i) * (2*pi) * EnrED_sp(i,j,k) / (lon) * (lon/ (2*pi)) );
            num2_edd_length(j)  = num2_edd_length(j)  + ( (EnrED_sp(i,j,k) * (lon/ (2*pi)) ) * (KX_INV(i)* ((2*pi) /lon)) ) ;
           
            
        end
    end
end

Eddy_length = (2*pi)*(num_edd_length ./ dnum_edd_length);
Eddy_length2 = (2*pi)*(num2_edd_length ./ num_edd_length);

%%%%%%%%%%%% END Transient eddy length scale


Enr_sp_av = zeros(1,bins);
Enr_spint_av = zeros(1,bins);
PEnr_sp_av = zeros(1,bins);
EnrED_sp_av = zeros(1,bins);
PEnrED_sp_av = zeros(1,bins);
MEnr_sp_av = zeros(1,bins);
MPEnr_sp_av = zeros(1,bins);

Enr_sp_av1 = zeros(1,bins);
PEnr_sp_av1 = zeros(1,bins);
EnrED_sp_av1 = zeros(1,bins);
PEnrED_sp_av1 = zeros(1,bins);
MEnr_sp_av1 = zeros(1,bins);
MPEnr_sp_av1 = zeros(1,bins);

Enr_sp_av2 = zeros(1,bins);
PEnr_sp_av2 = zeros(1,bins);
EnrED_sp_av2 = zeros(1,bins);
PEnrED_sp_av2 = zeros(1,bins);
MEnr_sp_av2 = zeros(1,bins);
MPEnr_sp_av2 = zeros(1,bins);

Enr_sp_av_LAY = zeros(bins,zsize);
Enr_sp_av1_LAY = zeros(bins,zsize);
Enr_sp_av2_LAY = zeros(bins,zsize);

EnrED_sp_av_LAY = zeros(bins,zsize);
EnrED_sp_av1_LAY = zeros(bins,zsize);
EnrED_sp_av2_LAY = zeros(bins,zsize);

sum_spc = 0;
sum_spc1 = 0;
sum_spc2 = 0;
for k=1:zsize
    for j=1:ysize
                if (Y(j) > -60) && (Y(j) < -40)
                    sum_spc = sum_spc +1;
                    Enr_sp_av(:) = Enr_sp_av(:) + (Enr_sp(:,j,k)); 
                    Enr_spint_av(:) = Enr_spint_av(i) + (Enr_spint(:,j,k));
                    PEnr_sp_av(:) = PEnr_sp_av(:) + ( PEnr_sp(:,j,k) * gprime(k));
                    EnrED_sp_av(:) = EnrED_sp_av(:) + ( EnrED_sp(:,j,k) );
                    PEnrED_sp_av(:) = PEnrED_sp_av(:) + ( PEnrED_sp(:,j,k) * gprime(k));
                    MEnr_sp_av(:) = MEnr_sp_av(:) + (MEnr_sp(:,j,k));
                    MPEnr_sp_av(:) = MPEnr_sp_av(:) + ( MPEnr_sp(:,j,k) * gprime(k));
                    Enr_sp_av_LAY(:,k) = Enr_sp_av_LAY(:,k) + (Enr_sp(:,j,k));
                    EnrED_sp_av_LAY(:,k) = EnrED_sp_av_LAY(:,k) + (EnrED_sp(:,j,k));
                end
                
                if (Y(j) >= -33.3)  %%% NORTH
                    sum_spc1 = sum_spc1 +1;
                    Enr_sp_av1(:) = Enr_sp_av1(:) + (Enr_sp(:,j,k) *(lon1 / lon));  
                    PEnr_sp_av1(:) = PEnr_sp_av1(:) + ( PEnr_sp(:,j,k) * gprime(k)   *(lon1 / lon));
                    EnrED_sp_av1(:) = EnrED_sp_av1(:) + ( EnrED_sp(:,j,k)  *(lon1 / lon) );
                    PEnrED_sp_av1(:) = PEnrED_sp_av1(:) + ( PEnrED_sp(:,j,k) * gprime(k)  *(lon1 / lon) );
                    MEnr_sp_av1(:) = MEnr_sp_av1(:) + (MEnr_sp(:,j,k) *(lon1 / lon));
                    MPEnr_sp_av1(:) = MPEnr_sp_av1(:) + ( MPEnr_sp(:,j,k) * gprime(k)  *(lon1 / lon) );
                    Enr_sp_av1_LAY(:,k) = Enr_sp_av1_LAY(:,k) + (Enr_sp(:,j,k) *(lon1 / lon));
                    EnrED_sp_av1_LAY(:,k) = EnrED_sp_av1_LAY(:,k) + (EnrED_sp(:,j,k) *(lon1 / lon));
            
                elseif (Y(j) < -33.3)  %%% SOUTH
                    sum_spc2 = sum_spc2 +1;
                    Enr_sp_av2(:) = Enr_sp_av2(:) + (Enr_sp(:,j,k)  *(lon2 / lon));
                    PEnr_sp_av2(:) = PEnr_sp_av2(:) + ( PEnr_sp(:,j,k) * gprime(k)  *(lon2 / lon));
                    EnrED_sp_av2(:) = EnrED_sp_av2(:) + ( EnrED_sp(:,j,k)   *(lon2 / lon) );
                    PEnrED_sp_av2(:) = PEnrED_sp_av2(:) + ( PEnrED_sp(:,j,k) * gprime(k)   *(lon2 / lon) );
                    MEnr_sp_av2(:) = MEnr_sp_av2(:) + (MEnr_sp(:,j,k)  *(lon2 / lon));
                    MPEnr_sp_av2(:) = MPEnr_sp_av2(:) + ( MPEnr_sp(:,j,k) * gprime(k)  *(lon2 / lon) );
                    Enr_sp_av2_LAY(:,k) = Enr_sp_av2_LAY(:,k) + (Enr_sp(:,j,k)  *(lon2 / lon));
                    EnrED_sp_av2_LAY(:,k) = EnrED_sp_av2_LAY(:,k) + (EnrED_sp(:,j,k)  *(lon2 / lon));
                end
           
        
    end
end

Enr_sp_av = Enr_sp_av /(sum_spc);
Enr_spint_av = Enr_spint_av /(sum_spc);
PEnr_sp_av = PEnr_sp_av /(sum_spc);
EnrED_sp_av = EnrED_sp_av /(sum_spc);
PEnrED_sp_av = PEnrED_sp_av /(sum_spc);
MEnr_sp_av = MEnr_sp_av /(sum_spc);
MPEnr_sp_av = MPEnr_sp_av /(sum_spc);


Enr_sp_av1 = Enr_sp_av1 /(sum_spc1);
PEnr_sp_av1 = PEnr_sp_av1 /(sum_spc1);
EnrED_sp_av1 = EnrED_sp_av1 /(sum_spc1);
PEnrED_sp_av1 = PEnrED_sp_av1 /(sum_spc1);
MEnr_sp_av1 = MEnr_sp_av1 /(sum_spc1);
MPEnr_sp_av1 = MPEnr_sp_av1 /(sum_spc1);

Enr_sp_av2 = Enr_sp_av2 /(sum_spc2);
PEnr_sp_av2 = PEnr_sp_av2 /(sum_spc2);
EnrED_sp_av2 = EnrED_sp_av2 /(sum_spc2);
PEnrED_sp_av2 = PEnrED_sp_av2 /(sum_spc2);
MEnr_sp_av2 = MEnr_sp_av2 /(sum_spc2);
MPEnr_sp_av2 = MPEnr_sp_av2 /(sum_spc2);

sum_spc_LAY = sum_spc / zsize;
sum_spc1_LAY = sum_spc1 / zsize;
sum_spc2_LAY = sum_spc2 / zsize;
Enr_sp_av_LAY = Enr_sp_av_LAY /(sum_spc_LAY);
Enr_sp_av1_LAY = Enr_sp_av1_LAY /(sum_spc1_LAY);
Enr_sp_av2_LAY = Enr_sp_av2_LAY /(sum_spc2_LAY);

EnrED_sp_av_LAY = EnrED_sp_av_LAY /(sum_spc_LAY);
EnrED_sp_av1_LAY = EnrED_sp_av1_LAY /(sum_spc1_LAY);
EnrED_sp_av2_LAY = EnrED_sp_av2_LAY /(sum_spc2_LAY);

num = num-1;
 
%%%%%%%%%%%% rms velocity
VEL_realta_GA1 = 0;
VEL_realta_GA2 = 0;
VEL_realta_GA1_LAY = zeros(zsize,1);
VEL_realta_GA2_LAY = zeros(zsize,1);
VEL_realtaED_GA1_LAY = zeros(zsize,1);
VEL_realtaED_GA2_LAY = zeros(zsize,1);
VEL_realta_GA1_ZY = zeros(ysize,zsize,1);
VEL_realta_GA2_ZY = zeros(ysize,zsize,1);
VEL_realtaED_GA1_ZY = zeros(ysize,zsize,1);
VEL_realtaED_GA2_ZY = zeros(ysize,zsize,1);

U_tta_za1 = zeros(ysize,1);
V_tta_za1 = zeros(ysize,1);
VEL_Tta_za1 = zeros(ysize,1);
VEL_Bta_za1 = zeros(ysize,1);

U_tta_za2 = zeros(ysize,1);
V_tta_za2 = zeros(ysize,1);
VEL_Tta_za2 = zeros(ysize,1);
VEL_Bta_za2 = zeros(ysize,1);


gsum_V1 = 0;
gsum_V2 = 0;
ARR_GA1 = 0;
ARR_GA2 = 0;

for i=1:xsize
    for j=1:ysize
        for k=1:zsize
                if (Y(j) >= -33.3)  %%% NORTH
                    gsum_V1 = gsum_V1 + 1;
                    VEL_realta_GA1 = VEL_realta_GA1 + VEL_realta(i,j,k) * (ARR(j));
                    VEL_realta_GA1_LAY(k) = VEL_realta_GA1_LAY(k) + VEL_realta(i,j,k) * (ARR(j));
                    VEL_realtaED_GA1_LAY(k) = VEL_realtaED_GA1_LAY(k) + VEL_realta_ED(i,j,k) * (ARR(j));    
                    VEL_realta_GA1_ZY(j,k) = VEL_realta_GA1_ZY(j,k) + VEL_realta(i,j,k) * (ARR(j));
                    VEL_realtaED_GA1_ZY(j,k) = VEL_realtaED_GA1_ZY (j,k) + VEL_realta_ED(i,j,k) * (ARR(j));
                    U_tta_za1(j)              = U_tta_za1(j) + U_tta(i,j) * ARR(j);
                    V_tta_za1(j)              = V_tta_za1(j) + V_tta(i,j) * ARR(j);
                    VEL_Tta_za1(j)              = VEL_Tta_za1(j) + VEL_Tta(i,j); % * ARR(j);  
                    VEL_Bta_za1(j)              = VEL_Bta_za1(j) + VEL_Bta(i,j); % * ARR(j);
                    ARR_GA1 = ARR_GA1 + ARR(j);
                elseif (Y(j) < -33.3)  %%% SOUTH
                    gsum_V2 = gsum_V2 + 1;
                    VEL_realta_GA2 = VEL_realta_GA2 + VEL_realta(i,j,k) * (ARR(j));
                    VEL_realta_GA2_LAY(k) = VEL_realta_GA2_LAY(k) + VEL_realta(i,j,k) * (ARR(j));
                    VEL_realtaED_GA2_LAY(k) = VEL_realtaED_GA2_LAY(k) + VEL_realta_ED(i,j,k) * (ARR(j));
                    VEL_realta_GA2_ZY(j,k) = VEL_realta_GA2_ZY(j,k) + VEL_realta(i,j,k) * (ARR(j));
                    VEL_realtaED_GA2_ZY(j,k) = VEL_realtaED_GA2_ZY (j,k) + VEL_realta_ED(i,j,k) * (ARR(j));
                    U_tta_za2(j)              = U_tta_za2(j) + U_tta(i,j) * ARR(j);
                    V_tta_za2(j)              = V_tta_za2(j) + V_tta(i,j) * ARR(j);
                    VEL_Tta_za2(j)              = VEL_Tta_za2(j) + VEL_Tta(i,j); %* ARR(j);
                    VEL_Bta_za2(j)              = VEL_Bta_za2(j) + VEL_Bta(i,j); % * ARR(j);
                    ARR_GA2 = ARR_GA2 + ARR(j);
                end 
        end
    end
end

ARR_GA = AR_GA;
gsum_V1_LAY = gsum_V1 / (zsize);
gsum_V2_LAY = gsum_V2 / (zsize);
ysize_N = gsum_V1 / (xsize * zsize); 
ysize_S = gsum_V2 / (xsize * zsize); 

gsum_V1_ZY = gsum_V1 / (ysize_N*zsize);
gsum_V2_ZY = gsum_V2 / (ysize_S*zsize);

gsum_V1_za = gsum_V1 / (ysize_N*zsize);
gsum_V2_za = gsum_V2 / (ysize_S*zsize);

%ARR_GA1 = ARR_GA1 / gsum_V1;
%ARR_GA2 = ARR_GA2 / gsum_V2;

VEL_realta_GA1 = (VEL_realta_GA1 /gsum_V1)/ ARR_GA;
VEL_realta_GA2 = (VEL_realta_GA2 /gsum_V2)/ ARR_GA;

VEL_realta_GA1_LAY = (VEL_realta_GA1_LAY /gsum_V1_LAY)/ ARR_GA;
VEL_realta_GA2_LAY = (VEL_realta_GA2_LAY /gsum_V2_LAY)/ ARR_GA;
VEL_realtaED_GA1_LAY = (VEL_realtaED_GA1_LAY /gsum_V1_LAY)/ ARR_GA;
VEL_realtaED_GA2_LAY = (VEL_realtaED_GA2_LAY /gsum_V2_LAY)/ ARR_GA;

VEL_realta_GA1_ZY = (VEL_realta_GA1_ZY /gsum_V1_ZY)/ ARR_GA;
VEL_realta_GA2_ZY = (VEL_realta_GA2_ZY /gsum_V2_ZY)/ ARR_GA;
VEL_realtaED_GA1_ZY = (VEL_realtaED_GA1_ZY /gsum_V1_ZY)/ ARR_GA;
VEL_realtaED_GA2_ZY = (VEL_realtaED_GA2_ZY /gsum_V2_ZY)/ ARR_GA;

U_tta_za1 = (U_tta_za1 / gsum_V1_za) / ARR_GA;
V_tta_za1 = (V_tta_za1 / gsum_V1_za) / ARR_GA;
VEL_Tta_za1 = (VEL_Tta_za1 / gsum_V1_za);
VEL_Bta_za1 = (VEL_Bta_za1 / gsum_V1_za);

U_tta_za2 = (U_tta_za2 / gsum_V2_za) / ARR_GA;
V_tta_za2 = (V_tta_za2 / gsum_V2_za) / ARR_GA;
VEL_Tta_za2 = (VEL_Tta_za2 / gsum_V2_za);
VEL_Bta_za2 = (VEL_Bta_za2 / gsum_V2_za);



dimen_num = 2* pi *(num / (90 * 111.32*1000*cosd(50)));
nd_dimen_num = 2* pi *(num / lon);
lombda = (2*pi) ./ (dimen_num);

dimen_num1 = 2* pi *(num / (90 * 111.32*1000*cosd(23.75)));
nd_dimen_num1 = 2* pi *(num / lon1);
dimen_num2 = 2* pi *(num / (90 * 111.32*1000*cosd(49.15)));
nd_dimen_num2 = 2* pi *(num / lon2);


lombda1 = (2*pi) ./ (dimen_num1);
lombda2 = (2*pi) ./ (dimen_num2);


length1 = size(lombda1,2);
length2 = size(lombda2,2);
V_rms1 = VEL_realta_GA1* ones(length1,1);
V_rms2 = VEL_realta_GA2 * ones(length2,1);

length1_LAY = ones(length1,zsize);
length2_LAY = ones(length2,zsize);
V_rms1_LAY  = ones(length2,zsize);
V_rms2_LAY  = ones(length2,zsize);

for k=1:zsize
    V_rms1_LAY(:,k) = VEL_realta_GA1_LAY(k) .* length1_LAY(:,k);
    V_rms2_LAY(:,k) = VEL_realta_GA2_LAY(k) .* length2_LAY(:,k);
end
    

VEL_edtaMR = zeros(ysize,1);
sum_VEL = 0;
for j=1:ysize
    for k=1:zsize
        for i=1:xsize
               sum_VEL = sum_VEL +1;
               VEL_edtaMR(j) = VEL_edtaMR(j)  + ( VEL_edta(i,j,k) );
        end
    end
end

sum_VEL = sum_VEL / ysize;
VEL_edtaMR = VEL_edtaMR / sum_VEL;
LR_ed = (VEL_edtaMR ./ beta).^0.5;


%%%%%%%%%%% SGS KE
ksize = size(Enr_sp_av,2);
KE_int =zeros(ksize,1);
DK = nd_dimen_num(10) - nd_dimen_num(9);

KE_tot = 0;
for i=1:ksize
    KE_tot = KE_tot + Enr_sp_av(i) * DK;
end

KE_sum =0;
for j=1:ksize
    KE_sum = KE_sum + Enr_sp_av(j) * DK;
    KE_same = Enr_sp_av(j) * DK;
    KE_int(j) = KE_tot - KE_sum + KE_same;
end
%%% NORTH
ksize1 = size(Enr_sp_av1,2);
KE_int1 =zeros(ksize1,1);
DK1 = nd_dimen_num1(10) - nd_dimen_num1(9);

KE_tot1 = 0;
for i=1:ksize1
    KE_tot1 = KE_tot1 + Enr_sp_av1(i) * DK1;
end

KE_sum1 =0;
for j=1:ksize1
    KE_sum1 = KE_sum1 + Enr_sp_av1(j) * DK1;
    KE_same1 = Enr_sp_av1(j) * DK1;
    KE_int1(j) = KE_tot1 - KE_sum1 + KE_same1;
end

%%%% SOUTH
ksize2 = size(Enr_sp_av2,2);
KE_int2 =zeros(ksize2,1);
DK2 = nd_dimen_num2(10) - nd_dimen_num2(9);

KE_tot2 = 0;
for i=1:ksize2
    KE_tot2 = KE_tot2 + Enr_sp_av2(i) * DK2;
end

KE_sum2 =0;
for j=1:ksize
    KE_sum2 = KE_sum2 + Enr_sp_av2(j) * DK2;
    KE_same2 = Enr_sp_av2(j) * DK2;
    KE_int2(j) = KE_tot2 - KE_sum2 + KE_same2;
end


%%%%%%%%%%%%%%%%%% SGS KE LAYER
ksize = size(Enr_sp_av_LAY,1);
KE_int_LAY =zeros(ksize,zsize);
KEED_int_LAY =zeros(ksize,zsize);
DK = nd_dimen_num(10) - nd_dimen_num(9);

Enr_sp_av_LAY(Enr_sp_av_LAY==0) = nan;
KE_tot = zeros(zsize,1);
KEED_tot = zeros(zsize,1);
for k=1:zsize
    for i=1:ksize
        if (~isnan(Enr_sp_av_LAY(i,k)))
            KE_tot(k) = KE_tot(k) + Enr_sp_av_LAY(i,k) * DK;
            KEED_tot(k) = KEED_tot(k) + EnrED_sp_av_LAY(i,k) * DK;
        end
    end
end

KE_sum =zeros(zsize,1);
KEED_sum =zeros(zsize,1);
KE_same =zeros(zsize,1);
KEED_same =zeros(zsize,1);
for k=1:zsize
    for j=1:ksize
        if (~isnan(Enr_sp_av_LAY(j,k)))
            KE_sum(k) = KE_sum(k) + Enr_sp_av_LAY(j,k) * DK;
            KEED_sum(k) = KEED_sum(k) + EnrED_sp_av_LAY(j,k) * DK;
            KE_same(k) = Enr_sp_av_LAY(j,k) * DK;
            KEED_same(k) = EnrED_sp_av_LAY(j,k) * DK;
            KE_int_LAY(j,k) = KE_tot(k) - KE_sum(k) + KE_same(k);
            KEED_int_LAY(j,k) = KEED_tot(k) - KEED_sum(k) + KEED_same(k);
        end
    end
end
Enr_sp_av_LAY(isnan(Enr_sp_av_LAY)) = 0;
%%% NORTH
ksize1 = size(Enr_sp_av1_LAY,1);
KE_int1_LAY =zeros(ksize1,zsize);
KEED_int1_LAY =zeros(ksize1,zsize);
DK1 = nd_dimen_num1(10) - nd_dimen_num1(9);

Enr_sp_av1_LAY(Enr_sp_av1_LAY==0) = nan;
KE_tot1 = zeros(zsize,1);
KEED_tot1 = zeros(zsize,1);
for k=1:zsize
    for i=1:ksize1
        if (~isnan(Enr_sp_av1_LAY(i,k)))
            KE_tot1(k) = KE_tot1(k) + Enr_sp_av1_LAY(i,k) * DK1;
            KEED_tot1(k) = KEED_tot1(k) + EnrED_sp_av1_LAY(i,k) * DK1;
        end
    end
end

KE_sum1 =zeros(zsize,1);
KEED_sum1 =zeros(zsize,1);
KE_same1 =zeros(zsize,1);
KEED_same1 =zeros(zsize,1);
for k=1:zsize
    for j=1:ksize1
        if (~isnan(Enr_sp_av1_LAY(j,k)))
            KE_sum1(k) = KE_sum1(k) + Enr_sp_av1_LAY(j,k) * DK1;
            KEED_sum1(k) = KEED_sum1(k) + EnrED_sp_av1_LAY(j,k) * DK1;
            KE_same1(k) = Enr_sp_av1_LAY(j,k) * DK1;
            KEED_same1(k) = EnrED_sp_av1_LAY(j,k) * DK1;
            KE_int1_LAY(j,k) = KE_tot1(k) - KE_sum1(k) + KE_same1(k);
            KEED_int1_LAY(j,k) = KEED_tot1(k) - KEED_sum1(k) + KEED_same1(k);
        end
    end
end
Enr_sp_av1_LAY(isnan(Enr_sp_av1_LAY)) = 0;
%%%% SOUTH
ksize2 = size(Enr_sp_av2_LAY,1);
KE_int2_LAY =zeros(ksize2,zsize);
KEED_int2_LAY =zeros(ksize2,zsize);
DK2 = nd_dimen_num2(10) - nd_dimen_num2(9);

Enr_sp_av2_LAY(Enr_sp_av2_LAY==0) = nan;
KE_tot2 = zeros(zsize,1);
KEED_tot2 = zeros(zsize,1);
for k=1:zsize
    for i=1:ksize2
        if (~isnan(Enr_sp_av2_LAY(i,k)))
            KE_tot2(k) = KE_tot2(k) + Enr_sp_av2_LAY(i,k) * DK2;
            KEED_tot2(k) = KEED_tot2(k) + EnrED_sp_av2_LAY(i,k) * DK2;
        end
    end
end

KE_sum2 =zeros(zsize,1);
KEED_sum2 =zeros(zsize,1);
KE_same2 =zeros(zsize,1);
KEED_same2 =zeros(zsize,1);
for k=1:zsize
    for j=1:ksize
        if (~isnan(Enr_sp_av2_LAY(j,k)))
            KE_sum2(k) = KE_sum2(k) + Enr_sp_av2_LAY(j,k) * DK2;
            KEED_sum2(k) = KEED_sum2(k) + EnrED_sp_av2_LAY(j,k) * DK2;
            KE_same2(k) = Enr_sp_av2_LAY(j,k) * DK2;
            KEED_same2(k) = EnrED_sp_av2_LAY(j,k) * DK2;
            KE_int2_LAY(j,k) = KE_tot2(k) - KE_sum2(k) + KE_same2(k);
            KEED_int2_LAY(j,k) = KEED_tot2(k) - KEED_sum2(k) + KEED_same2(k);
        end
    end
end
Enr_sp_av2_LAY(isnan(Enr_sp_av2_LAY)) = 0;

cd LAY

MatrixV(:,1) = dimen_num;
MatrixV(:,2) = Enr_sp_av;
dlmwrite('spc_FBC_7080.dat',MatrixV,'delimiter',' ','precision', 20);

MatrixV1(:,1) = dimen_num1;
MatrixV1(:,2) = Enr_sp_av1;
dlmwrite('spc_FBC1_7080.dat',MatrixV1,'delimiter',' ','precision', 20);

MatrixV2(:,1) = dimen_num2;
MatrixV2(:,2) = Enr_sp_av2;
dlmwrite('spc_FBC2_7080.dat',MatrixV2,'delimiter',' ','precision', 20);


V_rms1 = VEL_realta_GA1_LAY.^0.5;
V_rms2 = VEL_realta_GA2_LAY.^0.5;
V_rms1_ED = VEL_realtaED_GA1_LAY.^0.5;
V_rms2_ED = VEL_realtaED_GA2_LAY.^0.5;

V_rms1_SURF = V_rms1(1) * ones(size(lombda,2),1);
V_rms2_SURF = V_rms2(1) * ones(size(lombda,2),1);
V_rms1ED_SURF = V_rms1_ED(1) * ones(size(lombda,2),1);
V_rms2ED_SURF = V_rms2_ED(1) * ones(size(lombda,2),1);


dlmwrite('Vrms1_FBC_7080.dat',V_rms1,'delimiter',' ','precision', 20);
dlmwrite('Vrms2_FBC_7080.dat',V_rms2,'delimiter',' ','precision', 20);

dlmwrite('Vrms1ED_FBC_7080.dat',V_rms1_ED,'delimiter',' ','precision', 20);
dlmwrite('Vrms2ED_FBC_7080.dat',V_rms2_ED,'delimiter',' ','precision', 20);


CUM_VEL_SUR = KE_int_LAY(:,1).^0.5;
CUM_VEL1_SUR = KE_int1_LAY(:,1).^0.5;
CUM_VEL2_SUR = KE_int2_LAY(:,1).^0.5;

CUM_VELED_SUR = KEED_int_LAY(:,1).^0.5;
CUM_VELED1_SUR = KEED_int1_LAY(:,1).^0.5;
CUM_VELED2_SUR = KEED_int2_LAY(:,1).^0.5;

dlmwrite('lombda_FBC.dat',lombda,'delimiter',' ','precision', 20);
dlmwrite('lombda1_FBC.dat',lombda1,'delimiter',' ','precision', 20);
dlmwrite('lombda2_FBC.dat',lombda2,'delimiter',' ','precision', 20);

dlmwrite('CUV_FBC_7080.dat',CUM_VEL_SUR,'delimiter',' ','precision', 20);
dlmwrite('CUV1_FBC_7080.dat',CUM_VEL1_SUR,'delimiter',' ','precision', 20);
dlmwrite('CUV2_FBC_7080.dat',CUM_VEL2_SUR,'delimiter',' ','precision', 20);

dlmwrite('CUVED_FBC_7080.dat',CUM_VELED_SUR,'delimiter',' ','precision', 20);
dlmwrite('CUVED1_FBC_7080.dat',CUM_VELED1_SUR,'delimiter',' ','precision', 20);
dlmwrite('CUVED2_FBC_7080.dat',CUM_VELED2_SUR,'delimiter',' ','precision', 20);



cd ..


    
    


return

