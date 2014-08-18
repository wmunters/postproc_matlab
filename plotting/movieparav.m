
close all

for i=1:200

filename = strcat('field_inst_',num2str(i),'.csv');

u = load(filename); u = reshape(u(:,4),[80 192]);

xmesh = linspace(0,pi,192); zmesh = linspace(0,1,80);

xt = [0.2617 0.7851];    zt = [0.25 0.3];
xtow = xt+0.02;
zsin = .1 + .1*sin(2*pi*xmesh/pi);

figure(1)
hold on
plot(xmesh,zsin,'k','LineWidth',4);
plot([xt(1) xt(1)],[zt(1)-.05 zt(1)+.05],'k','LineWidth',2);
plot([xt(2) xt(2)],[zt(2)-.05 zt(2)+.05],'k','LineWidth',2);
plot([xtow(1) xtow(1)],[0.15 0.25],'k','LineWidth',2);
plot([xtow(2) xtow(2)],[0.2 0.3],'k','LineWidth',2);
pcolor(xmesh,zmesh,u); shading interp; axis equal; axis tight; caxis([0 10]); colorbar;
set(gca,'XTick',[]); set(gca,'YTick',[]);

saveas(gca,strcat(filename,'.png'));
set(gca,'FontSize',12);

end