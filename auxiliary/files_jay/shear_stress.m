function tau_xz = shear_stress(uw,Nx,Ny,Nz,zcc,hub_height);

uw_m(1:Nz) = 0.0;
for k = 1:Nz
for i = 1:Nx
  for j = 1:Ny
    uw_m(k) = uw_m(k)+uw(i,j,k);   
  end
end
end

 uw_m(:) = uw_m(:)/(Nx*Ny);
 
figure 
plot(-uw_m(:),zcc(:),'--black','LineWidth',1.5);
xlabel('$\it{\tau_{xz}/u_{*}^{2}}$','FontSize',20,'FontName','Times',...
       'Interpreter','latex')
ylabel('$\it{z/z_{h}}$','FontSize',20,'FontName','Times','Interpreter','latex');

 set(gca,'XTick',0.0:0.2:1.0,'FontSize',18,'FontName','Times')
 set(gca,'YTick',0:2:10,'FontSize',18,'FontName','Times')
 xlim([0.0 1.0])
 ylim([0 10])

line([0 1],[1-1/2 1-1/2],'LineStyle',':','Color','black','LineWidth',1.0)
line([0 1],[1+1/2 1+1/2],'LineStyle',':','Color','black','LineWidth',1.0)


end
