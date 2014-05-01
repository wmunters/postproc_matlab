close all
clear all

istart = 37.00;
istep  = 0.01;
iend   = 39.00;

turbines = false;


    [prefix,~] = strtok(path,';');
    prefix = strcat(prefix,'\');

Nx = 256;
Nx2 = 2*Nx;
Ny = Nx;
Ny2 = 2*Ny;
Nz = 80;
Nx32 = 1.5*Nx;
Ny32 = 1.5*Ny;

Lxfull = 2*pi;
Lx = Lxfull*(1-1/Nx);
Lyfull = Lxfull;
Ly = Lyfull*(1-1/Ny);

xmesh = linspace(0,Lx,Nx);
xmesh32 = linspace(0,Lx,Nx32);
xmeshbig = [xmesh xmesh+Lxfull];
ymesh = linspace(0,Ly,Ny);
ymesh32 = linspace(0,Ly,Ny32);
ymeshbig = [ymesh ymesh+Lyfull];
zmesh = linspace(1/Nz,1,Nz);
U = zeros(256,512,80);

xmin = 0; xmax = 2*Lx;
ymin = 0; ymax = 2*Ly;

dx= xmesh(2) - xmesh(1);
dy= ymesh(2) - ymesh(1);

if(turbines)
    WF = load('windfarm.setup');
end

% PLOTTING LOOP
%----------------------------------------------------------------------
for i=istart:istep:iend
    
    i
    timename = num2str(i);
    if(length(timename)==2) % geen cijfers na de komma
        timename = strcat(timename,'.0000');
    else
        while(length(timename)<7)
            timename = strcat(timename,'0');
        end
    end
       
    % Load fringe forces
    fu = load(strcat('force_u_t_',timename,'.dat')); fu = vect(fu); fu = reshape(fu,[Nx Ny]);
    fv = load(strcat('force_v_t_',timename,'.dat')); fv = vect(fv); fv = reshape(fv,[Nx Ny]);
    fw = load(strcat('force_w_t_',timename,'.dat')); fw = vect(fw); fw = reshape(fw,[Nx Ny]);
    
    % Load divergences
    div_vel = load(strcat('div_vel_t_',timename,'.dat')); div_vel = vect(div_vel); div_vel = reshape(div_vel,[Nx Ny]);
    div_rhs = load(strcat('div_rhs_t_',timename,'.dat')); div_rhs = vect(div_rhs); div_rhs = reshape(div_rhs,[Nx Ny]);
    div_force = load(strcat('div_force_t_',timename,'.dat')); div_force = vect(div_force); div_force = reshape(div_force,[Nx Ny]);
    
    figure(1)
    subplot(3,1,1); pcolor(xmesh,ymesh,fu'); fixfig;
    subplot(3,1,2); pcolor(xmesh,ymesh,fv'); fixfig;
    subplot(3,1,3); pcolor(xmesh,ymesh,fw'); fixfig;
    saveas(gcf,strcat(prefix,'forces_t_',timename,'.png'));
    
    figure(2) 
    pcolor(xmesh,ymesh,div_vel'); fixfig; caxis([-1 1]); title('Div velocity')
    saveas(gcf,strcat(prefix,'div_u_t_',timename,'.png'));
    
    figure(3)
    pcolor(xmesh,ymesh,div_rhs'); fixfig; title('Div rhs')
    saveas(gcf,strcat(prefix,'div_poisson_t_',timename,'.png'));
    
    figure(4)
    pcolor(xmesh,ymesh,div_force'); fixfig;title('Div force')
    saveas(gcf,strcat(prefix,'div_f_t_',timename,'.png'));
    
end

