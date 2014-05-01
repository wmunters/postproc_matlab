function V_mean = mean_farm_vel(u,v,Nx,Ny,Nz,zcc,z_h)

Um(1:Nz) = 0.0;
Vm(1:Nz) = 0.0;
for k = 1:Nz
for i = 1:Nx
  for j = 1:Ny
    Um(k) = Um(k)+u(i,j,k);   
    Vm(k) = Vm(k)+v(i,j,k);
  end
end
end

 Um(:) = Um(:)/(Nx*Ny);
 Vm(:) = Vm(:)/(Nx*Ny);  

figure 
semilogx(zcc(:),Um(:),'-black','LineWidth',1.5);
xlabel('$\it{z/z_{h}}$','FontSize',24,'FontName','Times','Interpreter','latex');
ylabel('$\it{\left\langle\overline{u}\right\rangle/u_{*,hi}}$','FontSize',20,'FontName','Times',...
       'Interpreter','latex');

 set(gca,'XTick',[0.1 1 10],'FontSize',18,'FontName','Times')
 set(gca,'YTick',4:2:16,'FontSize',18,'FontName','Times')
 set(gca,'FontSize',20,'FontName','Times')
 xlim([0.05 20])
 ylim([4 16])

line([1-1/2 1-1/2],[4 16],'LineStyle',':','Color','black','LineWidth',1.0)
line([1+1/2 1+1/2],[4 16],'LineStyle',':','Color','black','LineWidth',1.0)

% figure
% semilogx(zcc(:),Vm(:),'-black','LineWidth',1.5);
% 
% xlabel('$\it{z/z_{h}}$','FontSize',24,'FontName','Times','Interpreter','latex');
% ylabel('$\it{\left\langle\overline{v}\right\rangle/u_{*,hi}}$','FontSize',20,'FontName','Times',...
%        'Interpreter','latex');
% 
%  set(gca,'XTick',[0.1 1 10],'FontSize',18,'FontName','Times')
%  set(gca,'YTick',4:2:16,'FontSize',18,'FontName','Times')
%  set(gca,'FontSize',20,'FontName','Times')
%  xlim([0.05 20])
%  ylim([4 16])
% 
% line([1-1/2 1-1/2],[4 16],'LineStyle',':','Color','black','LineWidth',1.0)
% line([1+1/2 1+1/2],[4 16],'LineStyle',':','Color','black','LineWidth',1.0)


end







