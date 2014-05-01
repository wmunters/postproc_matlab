clc
clear all

% read domain parameter
fid = fopen('NS.setup');
param = textscan(fid,'%s %s %*[^\n]',6);
fclose(fid);

Nx = str2num(param{1}{4}); Ny = str2num(param{2}{4});
Lx = str2num(param{1}{6}); Ly = str2num(param{2}{6});

% read windfarm.setup
fid = fopen('windfarm.setup');
farm_dat= textscan(fid,'%f %f %f %f %*[^\n]','delimiter','');
fclose(fid);

nturb = farm_dat{1}(1);
C_t   = farm_dat{1}(2);

farm_x = farm_dat{1}(3:end);
farm_y = farm_dat{2}(3:end);
farm_z = farm_dat{3}(3:end);
rad = farm_dat{4}(3:end);

diameter   = 2*rad(1);
hub_height = farm_z(1);

%read mesh file in wall normal direction
ZMESH_WF = ['ZMESH_WF80'];
zmesh = load(ZMESH_WF);
Nz  = (zmesh(1)-1)/2;
zst = zmesh(2:2:end)/hub_height;
zcc = zmesh(3:2:end)/hub_height;

% some normalizations
farm_x = farm_x/hub_height;
farm_y = farm_y/hub_height;
farm_z = farm_z/hub_height;

Lx = Lx/hub_height; Ly = Ly/hub_height; 
dx=Lx/(Nx-1); dy=Ly/(Ny-1);

%read BL fielsstat.dat
A=importdata('BL_fieldstat.dat');

nn = size(A,1)*size(A,2);
A  = reshape(A',[nn 1]);

  u = A(1:Nx*Ny*Nz);
  v = A(Nx*Ny*Nz+1:2*Nx*Ny*Nz);
  w = A(2*Nx*Ny*Nz+1:3*Nx*Ny*Nz);
  uu =A(3*Nx*Ny*Nz+1:4*Nx*Ny*Nz);
  vv =A(4*Nx*Ny*Nz+1:5*Nx*Ny*Nz);
  ww =A(5*Nx*Ny*Nz+1:6*Nx*Ny*Nz);
  uv =A(6*Nx*Ny*Nz+1:7*Nx*Ny*Nz);
  uw =A(7*Nx*Ny*Nz+1:8*Nx*Ny*Nz);
  vw =A(8*Nx*Ny*Nz+1:9*Nx*Ny*Nz);

  u = reshape(u,[Nx Ny Nz]);
  v = reshape(v,[Nx Ny Nz]);
  w = reshape(w,[Nx Ny Nz]);
  uu = reshape(uu,[Nx Ny Nz]);
  vv = reshape(vv,[Nx Ny Nz]);
  ww = reshape(ww,[Nx Ny Nz]);
  uv = reshape(uv,[Nx Ny Nz]);
  uw = reshape(uw,[Nx Ny Nz]);
  vw = reshape(vw,[Nx Ny Nz]);
  
% mean streamwise and spanwise velocit of the farm
 mean_farm_vel(u,v,Nx,Ny,Nz,zcc,hub_height);

 % stress profile
 shear_stress(uw,Nx,Ny,Nz,zcc,hub_height);
 
% mean velocities at different streamwise location
% profiles at spanwise locations can be plotted similarly
 time_average_vel(u,v,Nx,Ny,Nz,zcc,hub_height,dx,dy);
  
% turbulent intensity at different streamwise location
% profiles at spanwise locations can be plotted similarly
 turbulent_intensity(sqrt(uu),v,Nx,Ny,Nz,zcc,hub_height,dx,dy);

% color map of a streamwise velocity field
  color_map(u,zcc,hub_height,diameter,Lx,Ly,farm_x,farm_y)









