close all


k = 8;

save = true;
recycleplaneswitch = true;
plotfield = true;
plotline = true;
savepath = true;
stats = true;

if(savepath)
    [prefix,~] = strtok(path,';');
    prefix = strcat(prefix,'\');
else
    prefix='';
end&

Nx = 128;
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
%%%%%

xmesh = linspace(0,Lx,Nx);
xmesh32 = linspace(0,Lx,Nx32);
xmeshbig = [xmesh xmesh+Lxfull];
ymesh = linspace(0,Ly,Ny);
ymesh32 = linspace(0,Ly,Ny32);
ymeshbig = [ymesh ymesh+Lyfull];
zmesh = linspace(1/Nz,1,Nz);
U = zeros(256,512,80);

xmin = 0; xmax = 2*Lx;
ymin = Ly/2; ymax = 2*Ly;

dx= xmesh(2) - xmesh(1);
dy= ymesh(2) - ymesh(1);

startx  = 1;
starty  = 1;
cutoffx = Nx/2;
cutoffy = Ny/2;

startxp = 1;
startyp = 1;
cutoffxp = Nx2;
cutoffyp = Ny2;


alpha = 30*pi/180;

offsetplanes = [0; Lyfull];
hinge_span   = [Lx ; Ly]+offsetplanes;
hinge_stream   = [Lx ; Ly]+offsetplanes;
counter      = 1;

Nplot = 5;
iplot = [1 Nx/8 Nx/4 Nx/2 Nx*5/8];
colors = ['k' ;'b' ;'m'; 'r'; 'c'];


% Rotation matrix
R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
iR = inv(R);

end_stream = [Lx*(1-sin(alpha)) ; Ly - Lx*(cos(alpha))]+offsetplanes;
end_span   = [Lx+Ly*cos(alpha)     ; Ly*( 1 - sin(alpha))]+offsetplanes;
end_corner = [end_stream(1) + (end_span(1) - hinge_span(1));end_stream(2) + (end_span(2) - hinge_span(2))];


for k=1:80

k    
    
uz = BLst.u(:,:,k);
vz = BLst.v(:,:,k);
Umain = uz*iR(1,1) + vz*iR(1,2);
Vmain = uz*iR(2,1) + vz*iR(2,2);
Wmain = BLst.w(:,:,k);
pmain = BLst.p(:,:,k);

pmain_dummy = [pmain;pmain]; pmain_dummy =[pmain_dummy pmain_dummy];
[dpdymain,dpdxmain] = gradient(pmain_dummy,dx,dy);
dpdxmain = dpdxmain(Nx+1:end,Ny+1:end);
dpdymain = dpdymain(Nx+1:end,Ny+1:end);
clear pmain_dummy

uzp = BLstp.u(:,:,k);
vzp = BLstp.v(:,:,k);
Uprec_big = [uzp; uzp]; Uprec_big = [Uprec_big Uprec_big];
Vprec_big = [vzp; vzp]; Vprec_big = [Vprec_big Vprec_big];
Wprec_big = [BLstp.w(:,:,k); BLstp.w(:,:,k)]; Wprec_big = [Wprec_big Wprec_big];
pprec_big = [BLstp.p(:,:,k); BLstp.p(:,:,k)]; pprec_big = [pprec_big pprec_big];

[dpdyprec_big,dpdxprec_big] = gradient(pprec_big,dx,dy);


figure(1) % U components
clf
hold on
    if(recycleplaneswitch)
        plot_recycleplanes(hinge_stream, end_stream, hinge_span, end_span, end_corner, cutoffx, cutoffy, Nx, Ny)
    end
    plot_prec_and_main(xmeshbig,startxp,cutoffxp,ymeshbig,startyp,cutoffyp,Uprec_big, xmesh,startx,cutoffx,ymesh,starty,cutoffy,Umain,end_stream, alpha)
    xlim([xmin xmax]); ylim([ymin ymax]);
    title('U')
    filename= strcat(prefix,'U_k',num2str(k),'_stat.png');
    if(save)
        saveas(gcf,filename);
    end


figure(2) % V component
clf
hold on
    if(recycleplaneswitch)
        plot_recycleplanes(hinge_stream, end_stream, hinge_span, end_span, end_corner, cutoffx, cutoffy, Nx, Ny)
    end
    plot_prec_and_main(xmeshbig,startxp,cutoffxp,ymeshbig,startyp,cutoffyp,Vprec_big, xmesh,startx,cutoffx,ymesh,starty,cutoffy,Vmain,end_stream, alpha)
    xlim([xmin xmax]); ylim([ymin ymax]);
    title('V')
    filename= strcat(prefix,'V_k',num2str(k),'_stat.png');
    if(save)
        saveas(gcf,filename);
    end
    
figure(3) % W component
clf
hold on    
    if(recycleplaneswitch)
        plot_recycleplanes(hinge_stream, end_stream, hinge_span, end_span, end_corner, cutoffx, cutoffy, Nx, Ny)
    end
    plot_prec_and_main(xmeshbig,startxp,cutoffxp,ymeshbig,startyp,cutoffyp,Wprec_big, xmesh,startx,cutoffx,ymesh,starty,cutoffy,Wmain,end_stream, alpha)
    xlim([xmin xmax]); ylim([ymin ymax]);
    title('W')
    filename= strcat(prefix,'W_k',num2str(k),'_stat.png');
    if(save)
        saveas(gcf,filename);
    end
    
figure(4) % pressure
clf
hold on
pprec_big = pprec_big - mean(mean(pprec_big));
pmain = pmain - mean(mean(pmain));
    if(recycleplaneswitch)
        plot_recycleplanes(hinge_stream, end_stream, hinge_span, end_span, end_corner, cutoffx, cutoffy, Nx, Ny)
    end
    plot_prec_and_main(xmeshbig,startxp,cutoffxp,ymeshbig,startyp,cutoffyp,pprec_big, xmesh,startx,cutoffx,ymesh,starty,cutoffy,pmain,end_stream, alpha)
    xlim([xmin xmax]); ylim([ymin ymax]);
    % Aanpassing maken aan de colorbar hier
    caxis([min(min(pprec_big)), max(max(pprec_big))])
    title('p')
    filename= strcat(prefix,'p_k',num2str(k),'_stat.png');
    if(save)
        saveas(gcf,filename);
    end
    
figure(5) % pressure gradient
clf
hold on
    if(recycleplaneswitch)
        plot_recycleplanes(hinge_stream, end_stream, hinge_span, end_span, end_corner, cutoffx, cutoffy, Nx, Ny)
    end
    plot_prec_and_main(xmeshbig,startxp,cutoffxp,ymeshbig,startyp,cutoffyp,dpdxprec_big, xmesh,startx,cutoffx,ymesh,starty,cutoffy,dpdxmain,end_stream, alpha)
    xlim([xmin xmax]); ylim([ymin ymax]);
    caxis([min(min(dpdxprec_big)), max(max(dpdxprec_big))])
    title('\partial p / \partial x')
    filename= strcat(prefix,'dpdx_k',num2str(k),'_stat.png');
    if(save)
        saveas(gcf,filename);
    end

figure(6) % pressure gradient
clf
hold on
    if(recycleplaneswitch)
        plot_recycleplanes(hinge_stream, end_stream, hinge_span, end_span, end_corner, cutoffx, cutoffy, Nx, Ny)
    end
    plot_prec_and_main(xmeshbig,startxp,cutoffxp,ymeshbig,startyp,cutoffyp,dpdyprec_big, xmesh,startx,cutoffx,ymesh,starty,cutoffy,dpdymain,end_stream, alpha)
    xlim([xmin xmax]); ylim([ymin ymax]);
    caxis([min(min(dpdyprec_big)), max(max(dpdyprec_big))])
    title('\partial p / \partial y')
    filename= strcat(prefix,'dpdy_k',num2str(k),'_stat.png');
    if(save)
        saveas(gcf,filename);
    end
    
end