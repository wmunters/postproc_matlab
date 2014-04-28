


% xmesh = linspace(0,2*pi*(1-1/256),256);
% ymesh = linspace(0,pi*(1-1/128),128);

figure
pcolor(xmesh,ymesh,BL.u(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('U')
saveas(gcf,'u.png');

figure
pcolor(xmesh,ymesh,BL.v(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('V')
saveas(gcf,'v.png');

figure
pcolor(xmesh,ymesh,BL.w(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('W')
saveas(gcf,'w.png');

figure
pcolor(xmesh,ymesh,BL.p(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('p')
saveas(gcf,'p.png');

figure
pcolor(xmesh,ymesh,BL.uu(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('uu')
saveas(gcf,'uu.png');

figure
pcolor(xmesh,ymesh,BL.vv(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('vv')
saveas(gcf,'vv.png');

figure
pcolor(xmesh,ymesh,BL.ww(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('ww')
saveas(gcf,'ww.png');

figure
pcolor(xmesh,ymesh,BL.uv(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('uv')
saveas(gcf,'uv.png');

figure
pcolor(xmesh,ymesh,BL.uw(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('uw')
saveas(gcf,'uw.png');

figure
pcolor(xmesh,ymesh,BL.vw(:,:,8)')
axis equal; axis tight; colorbar; shading interp;
title('vw')
saveas(gcf,'vw.png');


close all


