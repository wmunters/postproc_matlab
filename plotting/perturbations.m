% Script for generating perturbations which will be superposed on inflow
% profiles in order to assess turbulence development from pseudo-laminar
% inflow profiles. 
close all
clear all

Lx = 4*pi;  Nx = 512;   Lxprime = Lx*(1-1/Nx);  xmesh = linspace(0,Lxprime,Nx);  
Ly = pi;    Ny = 128;   Lyprime = Ly*(1-1/Ny);  ymesh = linspace(0,Lyprime,Ny);
Lz = 1;     Nz = 80;    Lzprime = Lz*(1-1/Nz);  zmesh = linspace(0,Lzprime,Nz);

U = zeros(Nx,Ny,Nz);
ampl = 2;

dx = 10; dy = dx;
xstart = Nx-dx;
xstop  = Nx;
ystart = Ny-dy;
ystop  = Ny;
zmax   = 0.8;

Px = 0.5;
Py = 0.25;
Pz = 0.25;

for k=1:Nz
    z = zmesh(k);
    ampl_z(k) = ampl*cos(pi/2*z/zmax);
    ampl_z(k) = max(ampl_z(k), 0);

    % spanwise boundary
    for i=1:Nx
        for j=ystart:ystop
        x = Lx - xmesh(i); y = Ly - jymesh(j); 
    
        U(i,j,k) = ampl_z(k)*cos(2*pi*x/Px)*cos(2*pi*y/Py)*cos(2*pi*z/Pz);
%        U(i,j,k) = ampl_z(k)*cos(2*pi*x/Px);
        end
    end

    %streamwise boundary
    for i=xstart:xstop
        for j=1:ystart
        x = xmesh(i); y = ymesh(j); 
    
        U(i,j,k) = ampl_z(k)*cos(2*pi*x/Px)*cos(2*pi*y/Py)*cos(2*pi*z/Pz);
%        U(i,j,k) = ampl_z(k)*cos(2*pi*y/Py);
        end
    end

end

figure
subplot(3,1,1); pcolor(ymesh,zmesh,squeeze(U(end,:,:))'); shading interp; title('Side view (streamwise boundary)'); 
subplot(3,1,2); pcolor(xmesh,zmesh,squeeze(U(:,end,:))'); shading interp; title('Side view (spanwise boundary)'); 
subplot(3,1,3); pcolor(xmesh,ymesh,squeeze(U(:,:,8))'); shading interp; title('Top view')
