close all
clear all

istart = 2.34;
istep  = 0.01;
iend   = 2.34;

startcompx = 75;
startcompy = 75;
stopcompx = 50;
stopcompy  = 100;

perc_fringex = 0.80;
perc_fringey = 0.80;

save = true;
recycleplaneswitch = false;
plotfringe=true;
fastrun = false;
plotfield = true;
plotline = false;
plotinflowstream = false;
savepath = true;
plotpressure = false;
plotfieldzerodeg = false;
plotvar='hor';
turbines = false;

if(savepath)
    [prefix,~] = strtok(path,';');
    prefix = strcat(prefix,'\');
else
    prefix='';
end

Nx = 256;
Nx2 = 2*Nx;
Ny = 256;
Ny2 = 2*Ny;
Nz = 80;

Nx32 = 1.5*Nx;
Ny32 = 1.5*Ny;

rescaling_factor = cos(30*pi/180);




Lxfull = pi;
Lx = Lxfull*(1-1/Nx);
Lyfull = pi;
Ly = Lyfull*(1-1/Ny);
Lxfringe = perc_fringex*Lxfull;
Lyfringe = perc_fringey*Lyfull;
%%%%%

xrec = 0.75*Lx;

xmesh = linspace(0,Lx,Nx);
xmesh32 = linspace(0,Lx,Nx32);
xmeshbig = [xmesh xmesh+Lxfull xmesh+2*Lxfull xmesh+3*Lxfull xmesh+4*Lxfull];
ymesh = linspace(0,Ly,Ny);
ymesh32 = linspace(0,Ly,Ny32);
ymeshbig = [ymesh ymesh+Lyfull ymesh+2*Lyfull ymesh+3*Lyfull ymesh+4*Lyfull];
zmesh = linspace(1/Nz,1,Nz);
U = zeros(256,512,80);

xmin = -0.25*Lx; xmax = 1.25*Lx;
ymin = -0.25*Ly; ymax = 1.25*Ly;


xmin = -0.45*Lx; xmax = 1.25*Lx;
ymin = -0.45*Ly; ymax = 1.25*Ly;

dx= xmesh(2) - xmesh(1);
dy= ymesh(2) - ymesh(1);

startx  = 1;
starty  = 1;
cutoffx = Nx;
cutoffy = Ny;

startxp = 1;
startyp = 1;
cutoffxp = 4*Nx;
cutoffyp = 4*Ny;

if(plotvar=='v')
    cmin=-5;
    cmax=5;
else
    cmin=15; %9
    cmax=26; %22.5
end

 cmin=13
 cmax=23

cminp = -10
cmaxp = 10
cmindp = -300
cmaxdp = 300


%   
theta = load('alpha.dat');
theta(:,2) = theta(:,2)*pi/180; % THETA IS IN RADIANS, ALPHA IS IN DEGS

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
    f1=figure('Visible','off')
else
    f1=figure
end
% theta = load('fringe.dat');

x0 = 0; 
y0 = Ly;
xstart = 2*Lxfull;
ystart = 3*Lyfull;
Lprime = sqrt(xstart^2 + ystart^2);
alpha0 = atan(ystart/xstart);
N = 4;
% PLOTTING LOOP
%----------------------------------------------------------------------
for i=istart:istep:iend
    
    i
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
%        alpha = 0;
    alpha =30 *pi/180;
    
    
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
    % Load precursor velocity data
    uzp = load(strcat('u_prec_zplane_k008_t_',timename,'.dat')); uzp = vect(uzp); uzp = reshape(uzp,[Nx Ny]);
    vzp = load(strcat('v_prec_zplane_k008_t_',timename,'.dat')); vzp = vect(vzp); vzp = reshape(vzp,[Nx Ny]);

    % Determine variable to be plotted
    if(plotvar=='u')
         plotvar_main = uz.*iR(1,1) + vz.*iR(1,2);
        plotvar_prec = uzp;
    elseif(plotvar=='v')
        plotvar_main = uz.*iR(2,1) + vz.*iR(2,2);
        plotvar_prec = vzp;
    elseif(plotvar=='hor')
        plotvar_main = sqrt(uz.^2 + vz.^2);
        plotvar_prec = sqrt(uzp.^2 + vzp.^2);
    else
        disp('Error in selection of plotvariable')
        return
    end
    
        
    % Copy precursor data to domain 4 times original size
    plotvar_prec_big = [plotvar_prec; plotvar_prec; plotvar_prec; plotvar_prec; plotvar_prec]; plotvar_prec_big = [plotvar_prec_big plotvar_prec_big plotvar_prec_big plotvar_prec_big plotvar_prec_big];
        
    phi(counter) = mean(mean( atan(vz(startcompx:stopcompx,startcompy:stopcompy)./uz(startcompx:stopcompx,startcompy:stopcompy))));
    
        clf
        if(plotfield)
        % Planes
        
        if(recycleplaneswitch)
            % Full main domain
            plot([hinge_stream(1) end_stream(1)],[hinge_stream(2) end_stream(2)],'k','LineWidth',0.3)
            plot([hinge_span(1) end_span(1)],[hinge_span(2) end_span(2)],'k','LineWidth',0.3)
            plot([end_corner(1) end_span(1)],[end_corner(2) end_span(2)],'k','LineWidth',0.3)
            plot([end_corner(1) end_stream(1)],[end_corner(2) end_stream(2)],'k','LineWidth',0.3)
            % Add boundaries of smaller domain
            if(cutoffx ~= Nx || cutoffy ~= Ny)
                plot([downright_small(1) end_stream(1)],[downright_small(2) end_stream(2)],':k','LineWidth',0.5)
                plot([topleft_small(1) end_stream(1)],[topleft_small(2) end_stream(2)],':k','LineWidth',0.5)
                plot([topleft_small(1) topright_small(1)],[topleft_small(2) topright_small(2)],':k','LineWidth',0.5)
                plot([downright_small(1) topright_small(1)],[downright_small(2) topright_small(2)],':k','LineWidth',0.5)
            end
        end
       
        hold on

        % Precursor
        g=pcolor(xmeshbig(startxp:cutoffxp),ymeshbig(startyp:cutoffyp),plotvar_prec_big(startxp:cutoffxp,startyp:cutoffyp)');
        set(g,'FaceColor','interp')
        
        
        
         
        
        caxis manual
        caxis([cmin cmax])
%         cbar = colorbar;
%         lab = ylabel(cbar,'(u^2 + v^2)^{1/2}/u_\tau');
        shading interp
        axis equal
        axis tight
        
         for j=1:N
         for i=1:N
          f=patch([0 0 Lx Lx]+(i-1)*[Lx Lx Lx Lx],[0 Ly Ly 0]+(j-1)*[Ly Ly Ly Ly],'w');set(f,'FaceAlpha',0,'EdgeColor','w','LineWidth',2.5);
          g = [g;f];
         end
         end

        % Main
        
        
        
         
%          shading interp
%         f2=patch([0 0 Lx Lx]+2*[Lx Lx Lx Lx],[0 Ly Ly 0],'w');set(f2,'FaceAlpha',1);
%         f = rectangle('position',[0 0 Lx Ly],'edgecolor','k');
%         g = [g;f];
        
        if(cutoffx~=Nx)
            plot([xmesh(startx) xmesh(startx)],[ymesh(starty) ymesh(cutoffy)],':k')
            plot([xmesh(cutoffx) xmesh(cutoffx)],[ymesh(starty) ymesh(cutoffy)],':k')
            plot([xmesh(startx) xmesh(cutoffx)],[ymesh(starty) ymesh(starty)],':k')
            plot([xmesh(startx) xmesh(cutoffx)],[ymesh(cutoffy) ymesh(cutoffy)],':k')
        end
                

    	% Rotate Precursor
        for i=1:numel(g)
        rotatefig(g(i),alpha,x0,y0,Lprime,alpha0)
        end
%           rotaterect(f,alpha,x0,y0,Lprime,alpha0)
%          rotatefig(f2,alpha,x0,y0,Lprime,alpha0)
        h=pcolor(xmesh(startx:cutoffx),ymesh(starty:cutoffy),plotvar_main(startx:cutoffx,starty:cutoffy)');
        set(h,'ZData',0.01*ones(256),'EdgeColor','None','FaceColor','interp');
        % Reposition Main
        xm = get(h,'XData');
        ym = get(h,'Ydata');
        % Vertical alignment
        shift_y = end_stream(2);
        shift_y = 0;
        ym = ym + shift_y;
        set(h,'YData',ym);
        % Horizontal alignment
        shift_x = end_stream(1);
        shift_x = 0;
        xm = xm + shift_x;
        set(h,'XData',xm);
        
        set(gca,'xtick',[],'ytick',[]);
        xlim([xmin xmax]); ylim([ymin ymax]);
        
        if(turbines)
                WF2 = WF;
                WF2(:,1) = WF2(:,1) + Lx;
                WF2(:,2) = WF2(:,2) + Ly;
                plot_turbines_topview(WF2)
        end
        
        
        
        % Calculate flow direction in a certain box in downstream region 
        
%         plot([xmesh(startcompx) xmesh(startcompx)],[ymesh(startcompy) ymesh(stopcompy)],'k')
%         plot([xmesh(stopcompx) xmesh(stopcompx)],[ymesh(startcompy) ymesh(stopcompy)],'k')
%         plot([xmesh(startcompx) xmesh(stopcompx)],[ymesh(startcompy) ymesh(startcompy)],'k')
%         plot([xmesh(startcompx) xmesh(stopcompx)],[ymesh(stopcompy) ymesh(stopcompy)],'k')

        a = 1.5;
%         comp = compass(a*cos(alpha), a*sin(alpha));
%         set(comp,'Color','k');
%         comp2 = compass(a*cos(phi(counter)), a*sin(phi(counter)));
%         set(comp,'LineWidth',2)
%         set(comp2,'LineWidth',2)
%        set(comp2,'Color','k')
% %         sx = (xmesh(startcompx)+xmesh(stopcompx))/2;
%         sy = (ymesh(startcompy)+ymesh(stopcompy))/2;
%         set(comp2,'XData',get(comp2,'XData')+[sx sx sx sx sx])
%         set(comp2,'YData',get(comp2,'YData')+[sy sy sy sy sy])
        
        

        % Plot fringe region outline
        plotfringe = false;
        if(plotfringe)
            plot([Lxfringe Lxfringe],[0 Lyfringe],'-.k','LineWidth',1.5)
            plot([0 Lxfringe],[Lyfringe Lyfringe],'-.k','LineWidth',1.5)
        end
        
        plotregions = false;
        lw = 2.5;
        linestyle = '-c';
        linestylesp = '-.g';
        linestylec = '-ow';
        fills = 'r';
        fillsp = 'm';
        fillc = [1 0.5 0.2];
        opac = 0.6;
        zh = 0.1;
        %Streamwise fringe region
        lw = 2.5;
        ls = '-.';
        s=fill3([Lxfringe Lx Lx Lxfringe],[0 0 Lyfringe Lyfringe],[zh zh zh zh],fills,'LineWidth',lw); set(s,'FaceAlpha',opac);
        s2=fill3([Lxfringe Lx Lx Lxfringe]-[Lx Lx Lx Lx],[0 0 Lyfringe Lyfringe],[zh zh zh zh],fills,'LineStyle',ls,'LineWidth',lw); set(s2,'FaceAlpha',opac);
        %Spanwise fringe region
        sp=fill3([0 Lxfringe Lxfringe 0],[Lyfringe Lyfringe Ly Ly],[zh zh zh zh],fillsp,'LineWidth',lw); set(sp,'FaceAlpha',opac);
        sp2=fill3([0 Lxfringe Lxfringe 0],[Lyfringe Lyfringe Ly Ly]-[Ly Ly Ly Ly],[zh zh zh zh],fillsp,'LineStyle',ls,'LineWidth',lw); set(sp2,'FaceAlpha',opac);
        
        %Corner region
        c=fill3([Lxfringe Lx Lx Lxfringe],[Lyfringe Lyfringe Ly Ly],[zh zh zh zh],fillc,'LineWidth',lw);set(c,'FaceAlpha',opac);
        c2=fill3([Lxfringe Lx Lx Lxfringe]-[Lx Lx Lx Lx],[Lyfringe Lyfringe Ly Ly]-[Ly Ly Ly Ly],[zh zh zh zh],fillc,'EdgeColor','None');set(c2,'FaceAlpha',opac);
        plot3([Lxfringe-Lx 0],[Lyfringe-Ly Lyfringe-Ly],[zh zh],'-.k','LineWidth',lw)
        plot3([Lxfringe-Lx Lxfringe-Lx],[Lyfringe-Ly 0],[zh zh],'-.k','LineWidth',lw)
        
        
        
        if(plotregions)
            %Streamwise fringe region
            plot([Lxfringe Lxfringe],[0 Lyfringe],linestyle,'LineWidth',lw)
            plot([Lxfringe Lx],[Lyfringe Lyfringe],linestyle,'LineWidth',lw)
            plot([Lx Lx],[0 Lyfringe],linestyle,'LineWidth',lw)
            plot([Lxfringe Lx],[0 0],linestyle,'LineWidth',lw)
            plot([Lxfringe-Lx Lxfringe-Lx],[0 Lyfringe],linestyle,'LineWidth',lw)
            plot([Lxfringe-Lx Lx-Lx],[Lyfringe Lyfringe],linestyle,'LineWidth',lw)
            plot([Lx-Lx Lx-Lx],[0 Lyfringe],linestyle,'LineWidth',lw)
            plot([Lxfringe-Lx Lx-Lx],[0 0],linestyle,'LineWidth',lw)
            %Spanwise fringe region
            plot([0 Lxfringe],[Lyfringe Lyfringe],linestylesp,'LineWidth',lw)
            plot([0 Lxfringe],[Ly Ly],linestylesp,'LineWidth',lw)
            plot([Lxfringe Lxfringe],[Lyfringe Ly],linestylesp,'LineWidth',lw)
            plot([0 0],[Lyfringe Ly],linestylesp,'LineWidth',lw)
            plot([0 Lxfringe],[Lyfringe Lyfringe]-[Ly Ly],linestylesp,'LineWidth',lw)
            plot([0 Lxfringe],[Ly Ly]-[Ly Ly],linestylesp,'LineWidth',lw)
            plot([Lxfringe Lxfringe],[Lyfringe Ly]-[Ly Ly],linestylesp,'LineWidth',lw)
            plot([0 0],[Lyfringe Ly]-[Ly Ly],linestylesp,'LineWidth',lw)
            %Corner region
            plot([Lxfringe Lxfringe],[Lyfringe Ly],linestylec,'LineWidth',lw)
            plot([Lx Lx],[Lyfringe Ly],linestylec,'LineWidth',lw)
            plot([Lxfringe Lx],[Ly Ly],linestylec,'LineWidth',lw)
            plot([Lxfringe Lx],[Lyfringe Lyfringe],linestylec,'LineWidth',lw)
            plot([Lxfringe Lxfringe]-[Lx Lx],[Lyfringe Ly]-[Ly Ly],linestylec,'LineWidth',lw)
            plot([Lx Lx]-[Lx Lx],[Lyfringe Ly]-[Ly Ly],linestylec,'LineWidth',lw)
            plot([Lxfringe Lx]-[Lx Lx],[Ly Ly]-[Ly Ly],linestylec,'LineWidth',lw)
            plot([Lxfringe Lx]-[Lx Lx],[Lyfringe Lyfringe]-[Ly Ly],linestylec,'LineWidth',lw)
        end
        
        % Plot main domain boundaries
        lw = 2.5;
        linestyle='b';
        plot3([0 0],[0 Ly],[1 1],linestyle,'LineWidth',lw)
        plot3([Lx Lx],[0 Ly],[1 1],linestyle,'LineWidth',lw)
        plot3([0 Lx],[0 0],[1 1],linestyle,'LineWidth',lw)
        plot3([0 Lx],[Ly Ly],[1 1],linestyle,'LineWidth',lw)
        
        % Finally, save the picture
        filename= strcat(prefix,plotvar,'_t_',timename,'.png');
        
        
        
        if(save)
            saveas(gcf,filename);
        end
        end
        
    counter = counter+1;
end
