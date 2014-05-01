function contour = color_map(u,zcc,z_h,diameter,Lx,Ly,farm_x,farm_y)

Nx = size(u,1);
Ny = size(u,2);
Nz = size(u,3);
nz_hub = find(zcc >= z_h/z_h,1);

%xy plane
uh = u(:,:,nz_hub);
% add extra row&column (periodic)d
  uh(:,Ny+1) = uh(:,1);
  uh(Nx+1,:) = uh(1,:); 

x = (0:1:Nx)/Nx * Lx;
y = (0:1:Ny)/Ny * Ly;
z = ones(length(y),length(x));

figure 
%surf(1:Nx,1:Ny,squeeze(u(1:Nx,1:Ny,nz_hub))')
surf(x,y,squeeze(uh)')
view(0,90)
shading interp
colorbar; 
xlabel('$\it{x/z_{h}}$','FontSize',20,'FontName','Times','Interpreter','latex');
ylabel('$\it{y/z_{h}}$','FontSize',20,'FontName','Times','Interpreter','latex');

  set(gca,'XTick',[0:10:60],'FontSize',15,'FontName','Times')
  set(gca,'YTick',[0:10:30],'FontSize',15,'FontName','Times')
  xlim([0 63])
  ylim([0 31.5])

  
end