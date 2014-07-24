close all
 clear all

istart = 3.00;
istep  = 0.01;
iend   = 39.00;

save = true;
recycleplaneswitch = true;
fastrun = false;
plotfield = false;
plotline = false;
plotinflowstream = false;
savepath = true;
plotpressure = true;
plotfieldzerodeg = true;
plotvar='u';
turbines = true;

if(savepath)
    [prefix,~] = strtok(path,';');
    prefix = strcat(prefix,'\');
else
    prefix='';
end

Nx = 256;
Nx2 = 2*Nx;
Ny = Nx/2;
Ny2 = 2*Ny;
Nz = 80;

Nx32 = 1.5*Nx;
Ny32 = 1.5*Ny;

rescaling_factor = cos(30*pi/180);





Lxfull = 2*pi;
Lx = Lxfull*(1-1/Nx);
Lyfull = Lxfull/2;
Ly = Lyfull*(1-1/Ny);
%%%%%

xrec = 0.75*Lx;
recycling_start=  0.6*Lx;
recycling_stop = 0.75*Lx;
fringe_start = 0.85*Lx;
fringe_stop = Lx;

xmesh = linspace(0,Lx,Nx);
xmesh32 = linspace(0,Lx,Nx32);
xmeshbig = [xmesh xmesh+Lxfull];
ymesh = linspace(0,Ly,Ny);
ymesh32 = linspace(0,Ly,Ny32);
ymeshbig = [ymesh ymesh+Lyfull];
zmesh = linspace(1/Nz,1,Nz);
U = zeros(256,512,80);

xmin = 0; xmax = 2*Lx;
ymin = 0; ymax = Ly;

dx= xmesh(2) - xmesh(1);
dy= ymesh(2) - ymesh(1);

startx  = 1;
starty  = 1;
cutoffx = Nx;
cutoffy = Ny;

startxp = 1;
startyp = 1;
cutoffxp = Nx2;
cutoffyp = Ny2;

if(plotvar=='v')
    cmin=-5;
    cmax=5;
else
    cmin=9;
    cmax=22.5;
end

cminp = -10
cmaxp = 10
cmindp = -300
cmaxdp = 300


%   
% theta = load('alpha.dat');
% theta(:,2) = theta(:,2)*pi/180; % THETA IS IN RADIANS, ALPHA IS IN DEGS

if(turbines)
    WF = load('windfarm.setup');
end


offsetplanes = [0; Lyfull];
% figure('Visible','off')
% hinge_stream = [Lx ; Ly];
% hinge_span   = [0 ; Ly];
hinge_span   = [Lx ; Ly]+offsetplanes;
hinge_stream   = [Lx ; Ly]+offsetplanes;
counter      = 1;

Nplot = 6;
iplot = [1 Nx/8 Nx/4 3*Nx/8 Nx/2 Nx*5/8];
colors = ['k' ;'b' ;'m'; 'r'; 'c'; 'g'];

% theta = load('alpha.dat');

if(fastrun)
    figure('Visible','off')
else
    figure
end
% theta = load('fringe.dat');

% PLOTTING LOOP
%----------------------------------------------------------------------
for i=istart:istep:iend
    
    i
    
    
    % 1. Initialization and loading variables
    % #######################################################
    
    % First correct timename so we can read in our files.
    timename = num2str(i);
    if(length(timename)==1) % geen cijfers na de komma
        timename = strcat(timename,'.0000');
    else
        while(length(timename)<6)
            timename = strcat(timename,'0');
        end
    end
    
    % Adjust data for recycleplanes so we can draw them and position main
%         alpha = theta(counter,2)%*pi/180;
%         alpha = 30*pi/180;
         alpha = 0;
%     alpha =30 *pi/180;
    %     alpha = min((i-33.0)/1,1)*30*pi/180;
    counter = counter+1;
    
    % Rotation matrix
    R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
    iR = inv(R);
    
    end_stream = [Lx*(1-sin(alpha)) ; Ly - Lx*(cos(alpha))]+offsetplanes;
    end_span   = [Lx+Ly*cos(alpha)     ; Ly*( 1 - sin(alpha))]+offsetplanes;
    end_corner = [end_stream(1) + (end_span(1) - hinge_span(1));end_stream(2) + (end_span(2) - hinge_span(2))];
    
    fracx = cutoffx/Nx;
    fracy = cutoffy/Ny;
    downright_small = end_stream*(1-fracx) + end_corner*fracx;
    topleft_small = hinge_stream*(fracy) + end_stream*(1-fracy);
    topright_small = topleft_small + downright_small - end_stream;
    
    % Load main velocity data
    uz = load(strcat('u_zplane_k008_t_',timename,'.dat')); uz = vect(uz); uz = reshape(uz,[Nx Ny]);
    vz = load(strcat('v_zplane_k008_t_',timename,'.dat')); vz = vect(vz); vz = reshape(vz,[Nx Ny]);
    
    % Determine variable to be plotted
    if(plotvar=='u')
%        plotvar_main = uz;
         plotvar_main = uz.*iR(1,1) + vz.*iR(1,2);
    elseif(plotvar=='v')
 %       plotvar_main = vz;
        plotvar_main = uz.*iR(2,1) + vz.*iR(2,2);
    elseif(plotvar=='hor')
        plotvar_main = sqrt(uz.^2 + vz.^2);
    else
        return
    end
    
    % Load precursor velocity data
    uzp = load(strcat('u_prec_zplane_k008_t_',timename,'.dat')); uzp = vect(uzp); uzp = reshape(uzp,[Nx Ny]);
    vzp = load(strcat('v_prec_zplane_k008_t_',timename,'.dat')); vzp = vect(vzp); vzp = reshape(vzp,[Nx Ny]);
    
    % Determine variable to be plotted
    if(plotvar=='u') %
%         plotvar_prec = uzp*cos(theta(counter,2)) + vzp*(-sin(theta(counter,2)));
        plotvar_prec = uzp;
    elseif(plotvar=='v')
%         plotvar_prec = uzp*sin(theta(counter,2)) + vzp*cos(theta(counter,2));
        plotvar_prec = vzp;
    elseif(plotvar=='hor')
        plotvar_prec = sqrt(uzp.^2 + vzp.^2);
    else
        return
    end
    
    % Copy precursor data to domain 4 times original size
    plotvar_prec_big = [plotvar_prec; plotvar_prec]; plotvar_prec_big = [plotvar_prec_big plotvar_prec_big];
   
      figure(1)
       clf
       subplot(2,1,2)
       hold on
       pcolor(xmesh,ymesh, plotvar_main')
       colorbar; caxis([cmin cmax]); shading interp; axis equal; axis tight;
       if(turbines)
                WF2 = WF;
                WF2(:,1) = WF2(:,1) + Lx;
                plot_turbines_topview(WF)
       end 
               plot([fringe_start fringe_start],[0 Ly],':k','LineWidth',2)     
       plot([fringe_stop fringe_stop],[0 Ly],':k','LineWidth',2)            
       
        subplot(2,1,1)
        hold on
        pcolor(xmesh,ymesh, plotvar_prec')
        colorbar; caxis([cmin cmax]); shading interp; axis equal; axis tight;
 
       plot([recycling_start recycling_start],[0 Ly],'-.k','LineWidth',2)
       plot([recycling_stop recycling_stop],[0 Ly],'-.k','LineWidth',2)
       plot([fringe_start fringe_start],[0 Ly],':k','LineWidth',2)     
       plot([fringe_stop fringe_stop],[0 Ly],':k','LineWidth',2)     
        
       
       filename = strcat(prefix,'field_merged_',plotvar,'_t_',timename,'.png');
        if(save)
            saveas(gcf,filename);
        end
        
        
        
    end
    
   