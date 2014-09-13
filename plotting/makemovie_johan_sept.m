close all
clear all

prec = false;

tstart = 40.00
tstep  = 0.01
tend   = 42.00

Nx = 768;
Ny = 768;
Nz = 192;

iplot = 501;
jplot = 534;
kplot = 16;

Lxfull = 10;
Lyfull = 10;
Lzfull = 1;

Lx = Lxfull*(1-1/Nx);
Ly = Lyfull*(1-1/Ny);
Lz = Lzfull;

x_fringe= 0.875*Lx;

xmesh = linspace(0,Lx,Nx);
xmeshbig = [xmesh xmesh+Lxfull];
ymesh = linspace(0,Ly,Ny);
zmesh = linspace(1/80,1,80);

U = zeros(Ny,Nx,Nz);Up = zeros(Ny,Nx,Nz);
Ubig = zeros(Ny,2*Nx,Nz);

startx  = 1;
starty  = 1;
cutoffx = Nx;
cutoffy = Ny;

[prefix,~] = strtok(path,';');
    prefix = strcat(prefix,'\');

% Load turbines
farm = load('windfarm.setup');
farm(:,1) = farm(:,1);
farm(:,2) = farm(:,2);
rad = 0.08;
hubh = 0.07;

cmin = 8;
cmax = 25;

figure
for t = tstart:tstep:tend

    % Display the time we're plotting at this moment
    t

    % First correct timename so we can read in our files.
    timename = num2str(t);
    if(length(timename)==2) % geen cijfers na de komma
        timename = strcat(timename,'.0000');
    else
        while(length(timename)<7)
            timename = strcat(timename,'0');
        end
    end
    
    iplotname = num2str(iplot);
    while(length(iplotname)<3)
        iplotname = strcat('0',iplotname);
    end
    jplotname = num2str(jplot);
    while(length(jplotname)<3)
        jplotname = strcat('0',jplotname);
    end
    kplotname = num2str(kplot);
    while(length(kplotname)<3)
        kplotname = strcat('0',kplotname);
    end

    % Load the main domain planes
    xname = strcat('u_xplane_i',iplotname,'_t_',timename,'.dat');
    ux = load(xname); ux = vect(ux); ux = reshape(ux,[Ny Nz]);
    yname = strcat('u_yplane_j',jplotname,'_t_',timename,'.dat');
    uy = load(yname); uy = vect(uy); uy = reshape(uy,[Nx Nz]);
    zname = strcat('u_zplane_k',kplotname,'_t_',timename,'.dat');
    uz = load(zname); uz = vect(uz); uz = reshape(uz,[Nx Ny]);

    % Load the precursor domain planes    
    if(prec)
    xnamep = strcat('u_prec_xplane_i',iplotname,'_t_',timename,'.dat');
    uxp = load(xnamep); uxp = vect(uxp); uxp = reshape(uxp,[Ny Nz]);
    ynamep = strcat('u_prec_yplane_j',jplotname,'_t_',timename,'.dat');
    uyp = load(ynamep); uyp = vect(uyp); uyp = reshape(uyp,[Nx Nz]);
    znamep = strcat('u_prec_zplane_k',kplotname,'_t_',timename,'.dat');
    uzp = load(znamep); uzp = vect(uzp); uzp = reshape(uzp,[Nx Ny]);
    end

    % Fill in the U array
    U(:,iplot,:) = ux;
    U(jplot,:,:) = uy;
    U(:,:,kplot)   = uz';
    
    % Do the actual plotting
    figure(1)
    clf
    hold on 
    slice(xmesh(startx:cutoffx),ymesh(starty:cutoffy),zmesh,U(starty:cutoffy,startx:cutoffx,:),[xmesh(iplot)],[ymesh(jplot)],[zmesh(kplot)]);
    for i=1:size(farm,1)
         plotCircle3D([farm(i,1) farm(i,2) farm(i,3)],[-1 0 0],0.05);
    end

    % Do some plot formatting
    xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')
    axis equal; axis tight; caxis manual; caxis ([cmin cmax]); shading interp
    hold off
    view(3)
    
    % Save to file
    filename= strcat(prefix,'u_t_',timename,'.png');
    saveas(gcf,filename);

end
