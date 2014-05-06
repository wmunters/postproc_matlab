close all
clear all

istart = 0.02;
istep  = 0.02;
iend   = 39.00;

perc_fringe = 0.8;

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
Ny = Nx;
Ny2 = 2*Ny;
Nz = 80;

Nx32 = 1.5*Nx;
Ny32 = 1.5*Ny;

rescaling_factor = cos(30*pi/180);




Lxfull = 2*pi;
Lx = Lxfull*(1-1/Nx);
Lyfull = Lxfull;
Ly = Lyfull*(1-1/Ny);
Lxfringe = perc_fringe*Lxfull;
Lyfringe = perc_fringe*Lyfull;
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
    cmin=9;
    cmax=22.5;
end

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
    figure('Visible','off')
else
    figure
end
% theta = load('fringe.dat');

x0 = 0; 
y0 = Ly;
xstart = 2*Lxfull;
ystart = 3*Lyfull;
Lprime = sqrt(xstart^2 + ystart^2);
alpha0 = atan(ystart/xstart);

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
        alpha = theta(counter,2)%*pi/180;
%        alpha = 0;
    %      alpha = 0;
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
    plotvar_prec_big = [plotvar_prec; plotvar_prec; plotvar_prec; plotvar_prec; plotvar_prec]; plotvar_prec_big = [plotvar_prec_big plotvar_prec_big plotvar_prec_big plotvar_prec_big plotvar_prec_big];
        
    
    % 2. Plot topviews of velocity field
    % #######################################################
    if(plotfield)
        f1=figure(1);
        clf
        
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
       
%         subplot(1,2,1)
        hold on
        
        % Precursor
        g=pcolor(xmeshbig(startxp:cutoffxp),ymeshbig(startyp:cutoffyp),plotvar_prec_big(startxp:cutoffxp,startyp:cutoffyp)')
%        set(g,'Visible','off')
%         contourf(xmeshbig(startxp:cutoffxp),ymeshbig(startyp:cutoffyp),plotvar_prec_big(startxp:cutoffxp,startyp:cutoffyp)')
        
        % plot additional lines
%        for k=1:Nplot
%             plot([xmeshbig(Nx+iplot(k)) xmeshbig(Nx+iplot(k))], [0 2*Ly], '-w', 'LineWidth', 2)
%        end
        
        % Main
        h=pcolor(xmesh(startx:cutoffx),ymesh(starty:cutoffy),plotvar_main(startx:cutoffx,starty:cutoffy)');
%         h=contourf(xmesh(startx:cutoffx),ymesh(starty:cutoffy),plotvar_main(startx:cutoffx,starty:cutoffy)');
%        set(h,'Visible','off') 
        caxis manual
        caxis([cmin cmax])
        cbar = colorbar;
%         set(get(cbar,'ylabel'),'String','u/u_\tau')
%         lab = ylabel(cbar,'\sqrt{u^2 + v^2}/u_\tau','Interpreter','Latex')
%         lab = ylabel(cbar,'$\sqrt{u^2 + v^2}$','interpreter','latex','FontSize',13)
        title('$\sqrt{u^2 + v^2}$','interpreter','latex','FontSize',14)
%         set(lab,'Rotation',0)
%         set(lab,'Position',get(lab,'Position')+[3 1 0]);
        shading interp
        axis equal
        axis tight
        % Rotate Main
%        rotate(h,[0 0 1],-alpha/pi*180,[0 0 0]);

        plot([0 0],[0 Ly],'k')
        plot([Lx Lx],[0 Ly],'k')
        plot([0 Lx],[0 0],'k')
        plot([0 Lx],[Ly Ly],'k')
        
        if(cutoffx~=Nx)
            plot([xmesh(startx) xmesh(startx)],[ymesh(starty) ymesh(cutoffy)],':k')
            plot([xmesh(cutoffx) xmesh(cutoffx)],[ymesh(starty) ymesh(cutoffy)],':k')
            plot([xmesh(startx) xmesh(cutoffx)],[ymesh(starty) ymesh(starty)],':k')
            plot([xmesh(startx) xmesh(cutoffx)],[ymesh(cutoffy) ymesh(cutoffy)],':k')
        end
                

    	% Rotate Precursor
        rotate(g,[0 0 1],alpha/pi*180,[0 0 0]);
        % Reposition Precursor
        xp = get(g,'XData');
        yp = get(g,'YData');
        % Vertical alignment
        shift_y = y0 - Lprime*sin(alpha0 + alpha);
        yp = yp + shift_y;
        set(g,'YData',yp);
        % Horizontal alignment
        shift_x = x0 - Lprime*cos(alpha0 + alpha);
        xp = xp + shift_x;
        set(g,'XData',xp);
        
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
        
        xlabel('x')
        ylabel('y')
        
        if(turbines)
                WF2 = WF;
                WF2(:,1) = WF2(:,1) + Lx;
                WF2(:,2) = WF2(:,2) + Ly;
                plot_turbines_topview(WF2)
            end
        
        xlim([xmin xmax]); ylim([ymin ymax]);
        
        
%       subplot(1,2,2)
        start = 1;
        stop  = 20;
        phi = mean(mean( atan(vz(start:stop,start:stop)./uz(start:stop,start:stop))));
%         plot([xmesh(start) xmesh(start)],[ymesh(start) ymesh(stop)],'k')
%         plot([xmesh(stop) xmesh(stop)],[ymesh(start) ymesh(stop)],'k')
%         plot([xmesh(start) xmesh(stop)],[ymesh(start) ymesh(start)],'k')
%         plot([xmesh(start) xmesh(stop)],[ymesh(stop) ymesh(stop)],'k')

        a = 1.5;
        comp = compass(a*cos(alpha), a*sin(alpha));
%         comp2 = compass(a*cos(phi), a*sin(phi));
        set(comp,'LineWidth',2)
%         set(comp2,'LineWidth',2)
%         set(comp2,'Color','k')
%         return
        
        % Plot fringe region outline
        if(plotfringe)
            plot([Lxfringe Lxfringe],[0 Lyfringe],':k')
            plot([0 Lxfringe],[Lyfringe Lyfringe],':k')
        end
        
        
        % Finally, save the picture
        filename= strcat(prefix,plotvar,'_t_',timename,'.png');
        
        if(save)
            saveas(gcf,filename);
        end
    end
    
        if(plotpressure)
        % Load pressure data
        p =  load(strcat('p_zplane_k008_t_',timename,'.dat'))     ; p = vect(p)  ; p =  reshape(p ,[Nx Ny]);
        pp = load(strcat('p_prec_zplane_k008_t_',timename,'.dat')); pp = vect(pp); pp = reshape(pp,[Nx Ny]);
        
        % Make bigger domains
        pbig = [p;p]; pbig = [pbig pbig];
        ppbig = [pp;pp]; ppbig = [ppbig ppbig];
        % Compute pressure gradients
        [dpdybig, dpdxbig] = gradient(pbig,dx,dy); dpdx = dpdxbig(Nx+1:end,Ny+1:end); dpdy = dpdybig(Nx+1:end,Ny+1:end);
        [dpdypbig, dpdxpbig] = gradient(ppbig,dx,dy); dpdxp = dpdxpbig(Nx+1:end,Ny+1:end); dpdyp = dpdypbig(Nx+1:end,Ny+1:end);
       
        figure(4)
        clf
        % Very ad-hoc implementation for aligned 0° grids
        dpdx_merged = [dpdxp; dpdx];
        dpdy_merged = [dpdyp; dpdy];
        p_merged = [pp; p];
        
%         subplot(3,1,1)
            hold on
            pcolor(xmeshbig,ymesh, p_merged')
            if(turbines)
                WF2 = WF;
                WF2(:,1) = WF2(:,1) + Lx;
                plot_turbines_topview(WF2)
            end
            title('Pressure')
            colorbar; caxis([cminp cmaxp]); shading interp; axis equal; axis tight;
        filename = strcat(prefix,'pressures_t_',timename,'.png');
        if(save)
            saveas(gcf,filename);
        end
            
%         subplot(3,1,2)
        figure(5)
             hold on
             pcolor(xmeshbig,ymesh, dpdx_merged')
             if(turbines)
                WF2 = WF;
                WF2(:,1) = WF2(:,1) + Lx;
                plot_turbines_topview(WF2)
            end
             title('\partial p / \partial x')
             colorbar; caxis([cmindp cmaxdp]); shading interp; axis equal; axis tight;
             
             filename = strcat(prefix,'dpdx_t_',timename,'.png');
        if(save)
            saveas(gcf,filename);
        end
%             
%         subplot(3,1,3)
        clf
             hold on
             pcolor(xmeshbig,ymesh, dpdy_merged')
             if(turbines)
                WF2 = WF;
                WF2(:,1) = WF2(:,1) + Lx;
                plot_turbines_topview(WF2)
            end
             title('\partial p / \partial y')
             colorbar; caxis([cmindp cmaxdp]); shading interp; axis equal; axis tight; 
             
        filename = strcat(prefix,'dpdy_t_',timename,'.png');
        if(save)
            saveas(gcf,filename);
        end
        
        
    end
    
    if(plotfieldzerodeg)
        
       figure(8)
       clf
       subplot(2,1,1)
       hold on
%        title('Streamwise velocity')
       plotvar_big = [plotvar_prec; plotvar_main];
       pcolor(xmeshbig,ymesh, plotvar_big')
       colorbar; caxis([cmin cmax]); shading interp; axis equal; axis tight;
       if(turbines)
                WF2 = WF;
                WF2(:,1) = WF2(:,1) + Lx;
                plot_turbines_topview(WF2)
            end 
       subplot(2,1,2)
       hold on
       pcolor(xmeshbig,ymesh, p_merged')
       if(turbines)
                WF2 = WF;
                WF2(:,1) = WF2(:,1) + Lx;
                plot_turbines_topview(WF2)
            end
%        title('Pressure')
       colorbar; caxis([cminp cmaxp]); shading interp; axis equal; axis tight;
        
        
       
       filename = strcat(prefix,'field_merged_',plotvar,'_t_',timename,'.png');
        if(save)
            I = getframe(f1)
            imwrite(I.cdata,filename,'PNG')
%             saveas(gcf,filename);
        end
        
        
        
    end

    
    
    % 3. Plot profiles
    % #######################################################
    if(plotline)
        
        % Planes perpendicular to x
        figure(2)
        clf
        hold on
        st = 1;
        
        % Load reference signal
%         [stream, span] = load_inflow(i,Nx,Ny); %%%%%%%%%%
%         inflow_stream_u = stream(:,:,1).*iR(1,1) + stream(:,:,2).*iR(1,2);
%         inflow_stream_v = stream(:,:,1).*iR(2,1) + stream(:,:,2).*iR(2,2);
%         inflow_span_u = span(:,:,1).*iR(1,1) + span(:,:,2).*iR(1,2);
%         inflow_span_v = span(:,:,1).*iR(2,1) + span(:,:,2).*iR(2,2);
%         
%         if(plotvar=='u')
%            plotvar_in_stream = inflow_stream_u;
%            plotvar_in_span = inflow_span_u;
%         elseif(plotvar=='v')
%            plotvar_in_stream = inflow_stream_v;
%            plotvar_in_span = inflow_span_v;
%         end
        
        
%          streamplane = load('recycleplane_stream_t_35.8000.dat');
%          streamplane = vect(streamplane); s(1:Ny32,1) = streamplane(1:Ny32); s(1:Ny32,2) = streamplane(Ny32+1:end)+Ly;
        s_stream(:,1) = linspace(end_stream(1), hinge_stream(1), Ny)+dx;
        s_stream(:,2) = linspace(end_stream(2), hinge_stream(2), Ny);
        shiftvector = [cos(alpha); -sin(alpha)];
        
        
        for k=1:Nplot
%         for k=1:1
%            xmesh_tilt = linspace(end_stream(1),hinge_stream(1),Ny32) + (iplot(k) - 1)*dx*shiftvector(1);
%             ymesh_tilt = linspace(end_stream(2),hinge_stream(2),Ny32) + (iplot(k) - 1)*dy*shiftvector(2);
            xmesh_tilt = s_stream(:,1) + (iplot(k) - 1)*dx*shiftvector(1);
            ymesh_tilt = s_stream(:,2) + (iplot(k) - 1)*dy*shiftvector(2);
        
            for j=1:Ny
                u_int(j) = bilinear_interpolation(plotvar_prec_big, xmeshbig, ymeshbig, xmesh_tilt(j), ymesh_tilt(j)) ;
            end
            
                        %     plot(uz(:,iplot(k)),colors(k),'LineWidth',2)
            %     plot(uzp(:,iplot(k)),strcat(':',colors(k)),'LineWidth',2)
            plot(ymesh,plotvar_main(iplot(k),:)+(k-1)*10,colors(k),'LineWidth',1) %%+1!!!!!!
%             plot(uzp(iplot(k),:)+(k-1)*10,ymesh,strcat(':',colors(k)),'LineWidth',1)
            plot(ymesh,u_int+(k-1)*10,strcat(':',colors(k)),'LineWidth',1)
            if(iplot(k)==1)
%                 plot(ymesh32(1:st:end),plotvar_in_stream(1:st:end,8),'b')
            end
        end
        xlabel('y')
        ylabel('U')
        title('Streamwise')
        xlim([0 Ly])
        
        filename= strcat(prefix,'xprofiles_',plotvar,'_t_',timename,'.png');
        if(save)
            saveas(gcf,filename);
        end
        
        
        % Planes perpendicular to y
        figure(3)
        clf
        hold on
        
%          spanplane = load('recycleplane_span_t_35.8000.dat');
%          spanplane = vect(spanplane); s(1:Nx32,1) = spanplane(1:Nx32); s(1:Nx32,2) = spanplane(Nx32+1:end)+Ly;
        s_span(:,1) = linspace(end_stream(1), end_corner(1), Nx);
        s_span(:,2) = linspace(end_stream(2), end_corner(2), Nx);
        shiftvector = [sin(alpha); cos(alpha)];
        
        for k=1:Nplot
            xmesh_tilt = s_span(:,1) + (iplot(k) - 1)*dx*shiftvector(1);
            ymesh_tilt = s_span(:,2) + (iplot(k) - 1)*dy*shiftvector(2);
%         
%             figure(1)
%             hold on
%             plot(xmesh_tilt,ymesh_tilt,'ow')
            
            for j=1:Ny
                u_int(j) = bilinear_interpolation(plotvar_prec_big, xmeshbig, ymeshbig, xmesh_tilt(j), ymesh_tilt(j)) ;
            end
             
            plot(xmesh,plotvar_main(:,iplot(k))+(k-1)*10,colors(k),'LineWidth',1)
            plot(xmesh-dx,u_int+(k-1)*10,strcat(':',colors(k)),'LineWidth',1)
            if(iplot(k)==1)
%                 plot(ymesh32(1:st:end),plotvar_in_span(1:st:end,8),'b')
            end
            
        end
        
        xlabel('x')
        ylabel('U')
        title('Spanwise')
        xlim([0 Lx])
        
        filename= strcat(prefix,'yprofiles_',plotvar,'_t_',timename,'.png');
        if(save)
            saveas(gcf,filename);
        end
        
%         pause
        
    end
    
    
    % 4. Plot inflow along with reference signal
    % #######################################################
    if(plotinflowstream)
        figure(4)
        clf
        hold on
        [stream, span] = load_inflow(i,Nx,Ny);
        
%         streamplane = load('recycleplane_stream_t_35.8000.dat');
%         streamplane = vect(streamplane); s(1:Ny32,1) = streamplane(1:Ny32); s(1:Ny32,2) = streamplane(Ny32+1:end);
        
        xmesh_tilt = linspace(end_stream(1),hinge_stream(1),Ny);
        ymesh_tilt = linspace(end_stream(2),hinge_stream(2),Ny);
        
        for j=1:Ny
           u_int(j) = bilinear_interpolation(plotvar_prec_big, xmeshbig, ymeshbig, xmesh_tilt(j), ymesh_tilt(j)) ;
        end
        
        
        plot(ymesh32,stream(:,8,1),'-.m')
%         plot(ymesh,uzp(256,:),'r')
        plot(ymesh,plotvar_main(end,:),'k')
        plot(ymesh,plotvar_prec(end,:),':r')
        
        title('Inflow velocity at x=0')
        ylabel('y')
        hl = legend('Reference','Main','Precursor');
        set(hl,'Location','South');
        set(hl,'Orientation','Horizontal');
%         pause
        filename = strcat(prefix,'inflow_stream_',plotvar,'_t_',timename,'.png');
        if(save)
            saveas(gcf,filename);
        end
        
%         figure(4)
%         clf
%         hold on
%         plot(xmesh32,span(:,8,1),'-.m')
%         plot(xmesh,uz(:,1),'k')
%         title('Span inflow at y=0')
%         
%         
%         filename = strcat('inflow_span_t_',timename,'.png');
%         if(save)
%             saveas(gcf,filename);
%         end
    end
    
%      pause
    
end
    
%     pause
    
    
