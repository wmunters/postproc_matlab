close all
clear all

tstart = 4.00
tstep  = 0.005
tend   = 6.00

Nx = 256;
Ny = 128;
Nz = 80;

iplot = 256
iplot1 = 178
jplot = 118
kplot = 8

Lxfull = 2*pi;
Lyfull = pi;
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
rad = 0.05;


min = 8;
max = 25;

% uz = load('u_zplane_k008_t_47.0000000.dat');
% uz = reshape(uz',[size(uz,1)*size(uz,2) 1]);
% uz = reshape(uz,[512 256]);
% 
% figure
% hold on
% pcolor(xmesh,ymesh,uz')
% for i=1:48
%     plot([farm(i,1) farm(i,1)],[farm(i,2)-rad farm(i,2)+rad],'k','LineWidth',3)
% end
% plot([0.4*Lx 0.4*Lx],[0 Ly],'-.k','LineWidth',2)
% plot([0.95*Lx 0.95*Lx],[0 Ly],'-k','LineWidth',2)
% axis equal
% axis tight
% 
% shading interp
% xlabel('x [km]')
% ylabel('y [km]')
% figure('units','normalized','outerposition',[0 0 1 1])
figure
for t = tstart:tstep:tend



t

clf
% for i=46.00:0.005:49.00
% Load data


% First correct timename so we can read in our files.
    timename = num2str(t);
    if(length(timename)==1) % geen cijfers na de komma
        timename = strcat(timename,'.0000');
    else
        while(length(timename)<6)
            timename = strcat(timename,'0');
        end
    end
    
    iplotname1 = num2str(iplot1);
    while(length(iplotname1)<3)
        iplotname1 = strcat('0',iplotname1);
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

    
    
xname = strcat('u_xplane_i',iplotname,'_t_',timename,'.dat');
ux = load(xname); ux = vect(ux); ux = reshape(ux,[Ny Nz]);

xname = strcat('u_xplane_i',iplotname1,'_t_',timename,'.dat');
ux1 = load(xname); ux1 = vect(ux1); ux1 = reshape(ux1,[Ny Nz]);

xnamep = strcat('u_prec_xplane_i',iplotname,'_t_',timename,'.dat');
uxp = load(xnamep); uxp = vect(uxp); uxp = reshape(uxp,[Ny Nz]);


yname = strcat('u_yplane_j',jplotname,'_t_',timename,'.dat');
uy = load(yname); uy = vect(uy); uy = reshape(uy,[Nx Nz]);
ynamep = strcat('u_prec_yplane_j',jplotname,'_t_',timename,'.dat');
uyp = load(ynamep); uyp = vect(uyp); uyp = reshape(uyp,[Nx Nz]);

zname = strcat('u_zplane_k',kplotname,'_t_',timename,'.dat');
uz = load(zname); uz = vect(uz); uz = reshape(uz,[Nx Ny]);
znamep = strcat('u_prec_zplane_k',kplotname,'_t_',timename,'.dat');
uzp = load(znamep); uzp = vect(uzp); uzp = reshape(uzp,[Nx Ny]);

U(:,iplot1,:) = ux1;
U(:,iplot,:) = ux;
U(jplot,:,:) = uy;
U(:,:,kplot)   = uz';

Up(:,iplot,:) = uxp;
Up(jplot,:,:) = uyp;
Up(:,:,kplot)   = uzp';

Ubig(:,1:Nx,:) = Up;
Ubig(:,Nx+1:2*Nx,:) = U;

% figure(1)
% clf
% hold on 
%  
% start = Nx-40;
% cut =  Nx+40;
% % slice(xmeshbig,ymesh,zmesh,Ubig,[xmeshbig(end)],[ymesh(jplot)],[zmesh(kplot)]);
% slice(xmeshbig(start:cut),ymesh,zmesh,Ubig(:,start:cut,:),[],[ymesh(jplot)],[zmesh(kplot)]);
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
% shading interp
% %  for i=1:size(farm,1)
%      for i=1:6
%      plotCircle3D([farm(i,1)+Lxfull farm(i,2) farm(i,3)],[-1 0 0],0.05);
%  end
% % plot3([Lxfull Lxfull],[0 0],[0.1 0.1],'ok','MarkerFaceColor','k')
% % plot3([Lxfull Lxfull],[ymesh(jplot-2) ymesh(jplot-2)],[0.1 0.1],'ok','MarkerFaceColor','k')
% plot3([Lxfull+.05 Lxfull+.05],[0 Ly],[0.1 0.1],':k','MarkerFaceColor','k')
% plot3([Lxfull+.05 Lxfull+.05],[ymesh(jplot) ymesh(jplot)],[0 1],':k','MarkerFaceColor','k')
% 
% plot3([Lxfull-.05 Lxfull-.05],[0 Ly],[0.1 0.1],':k','MarkerFaceColor','k')
% plot3([Lxfull-.05 Lxfull-.05],[ymesh(jplot) ymesh(jplot)],[0 1],':k','MarkerFaceColor','k')
% axis equal 
% axis tight
% caxis manual
% caxis ([min max])
% % h2 = colorbar('location','southoutside');
% % set(h2,'Position',[0.65 0.3 0.28 0.0347]);
% hold off
% view(3)
% % xlim([5 7])
% filename= strcat(prefix,'uzoom_t_',timename,'.png');
% saveas(gcf,filename);
% 
% 
% 
% 
% 
% end
% return

margin = 0.05;
figure(1)
clf
subplot_tight(2,1,2,margin)
hold on 
% slice(xmesh(startx:cutoffx),ymesh(starty:cutoffy),zmesh,U(starty:cutoffy,startx:cutoffx,:),[ xmesh(490)],[ymesh(178)],[zmesh(8)]);
 plot3([x_fringe x_fringe],[0 Ly],[zmesh(kplot) zmesh(kplot)],'-.k','LineWidth',2)
 
slice(xmesh(startx:cutoffx),ymesh(starty:cutoffy),zmesh,U(starty:cutoffy,startx:cutoffx,:),[xmesh(iplot1) xmesh(iplot)],[ymesh(jplot)],[zmesh(kplot)]);
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
shading interp
 for i=1:size(farm,1)
     plotCircle3D([farm(i,1) farm(i,2) farm(i,3)],[-1 0 0],0.05);
 end

axis equal 
axis tight
caxis manual
caxis ([min max])
% h2 = colorbar('location','southoutside');
% set(h2,'Position',[0.65 0.3 0.28 0.0347]);
hold off
view(3)
% filename= strcat(prefix,'u_t_',timename,'.png');
% saveas(gcf,filename);


% figure(2)
% clf
subplot_tight(2,1,1,margin)
hold on 
% slice(xmesh(startx:cutoffx),ymesh(starty:cutoffy),zmesh,U(starty:cutoffy,startx:cutoffx,:),[ xmesh(490)],[ymesh(178)],[zmesh(8)]);
slice(xmesh(startx:cutoffx),ymesh(starty:cutoffy),zmesh,Up(starty:cutoffy,startx:cutoffx,:),[ xmesh(iplot)],[ymesh(jplot)],[zmesh(kplot)]);
plot3([x_fringe x_fringe],[0 Ly],[zmesh(kplot) zmesh(kplot)],'-.k','LineWidth',2)
 
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
shading interp
% for i=1:36
%     plotCircle3D([farm(i,1) farm(i,2) farm(i,3)],[-1 0 0],0.05);
% end
axis equal 
axis tight
caxis manual
caxis ([min max])
h2 = colorbar('location','eastoutside');
% set(h2,'Position',[0.65 0.3 0.28 0.0347]);
set(h2,'Position',[0.85 0.22 0.03 0.6]);
hold off
view(3)


filename= strcat(prefix,'ujoin_t_',timename,'.png');
saveas(gcf,filename);

end




% end

