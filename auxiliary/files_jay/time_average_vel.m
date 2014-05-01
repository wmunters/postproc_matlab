

function V_time = time_average_vel(u,v,Nx,Ny,Nz,zcc,z_h,dx,dy)

% at the Ny/2 in spanwise direction
% 6 Nx position
nx_loc(1) = 1; nx_loc(2) = 25; nx_loc(3) = 50;
nx_loc(4) = 75; nx_loc(5) = 100; nx_loc(6) = 120;
nx_loc = nx_loc';

for i = 1:size(nx_loc,1);
  x_loc(i) = (nx_loc(i)-1)*dx;  
 for k = 1:Nz;
  u1(i,k) = u(nx_loc(i),Ny/2,k);
 end
end 
 
 
figure
 subplot(2,3,1)
 plot(u1(1,:),zcc(:),'-black','LineWidth',1.0);
 xlabel('$\it{u}/u_{*,hi}$','FontSize',16,'FontName','Times',...
       'Interpreter','latex');
 ylabel('$\it{z/z_{h}}$','FontSize',16,'FontName','Times','Interpreter','latex');


  set(gca,'XTick',[5:5:15],'FontSize',12,'FontName','Times')
  set(gca,'YTick',0:2:10,'FontSize',12,'FontName','Times')
  %set(gca,'FontSize',20,'FontName','Times')
  xlim([4 15])
  ylim([0 10])
 
 
 subplot(2,3,2)
 plot(u1(2,:),zcc(:),'-black','LineWidth',1.0);
 xlabel('$\it{u}/u_{*,hi}$','FontSize',16,'FontName','Times',...
       'Interpreter','latex');
 ylabel('$\it{z/z_{h}}$','FontSize',16,'FontName','Times','Interpreter','latex');


  set(gca,'XTick',[5:5:15],'FontSize',12,'FontName','Times')
  set(gca,'YTick',0:2:10,'FontSize',12,'FontName','Times')
  %set(gca,'FontSize',20,'FontName','Times')
  xlim([4 15])
  ylim([0 10])

 
 subplot(2,3,3)
 plot(u1(3,:),zcc(:),'-black','LineWidth',1.0);
 xlabel('$\it{u}/u_{*,hi}$','FontSize',16,'FontName','Times',...
       'Interpreter','latex');
 ylabel('$\it{z/z_{h}}$','FontSize',16,'FontName','Times','Interpreter','latex');


  set(gca,'XTick',[5:5:15],'FontSize',12,'FontName','Times')
  set(gca,'YTick',0:2:10,'FontSize',12,'FontName','Times')
  %set(gca,'FontSize',20,'FontName','Times')
  xlim([4 15])
  ylim([0 10]) 
 
 
 subplot(2,3,4) 
 plot(u1(4,:),zcc(:),'-black','LineWidth',1.0);
 xlabel('$\it{u}/u_{*,hi}$','FontSize',16,'FontName','Times',...
       'Interpreter','latex');
 ylabel('$\it{z/z_{h}}$','FontSize',16,'FontName','Times','Interpreter','latex');


  set(gca,'XTick',[5:5:15],'FontSize',12,'FontName','Times')
  set(gca,'YTick',0:2:10,'FontSize',12,'FontName','Times')
  %set(gca,'FontSize',20,'FontName','Times')
  xlim([4 15])
  ylim([0 10]) 
 
 
 subplot(2,3,5)
 plot(u1(5,:),zcc(:),'-black','LineWidth',1.0);
 xlabel('$\it{u}/u_{*,hi}$','FontSize',16,'FontName','Times',...
       'Interpreter','latex');
 ylabel('$\it{z/z_{h}}$','FontSize',16,'FontName','Times','Interpreter','latex');


  set(gca,'XTick',[5:5:15],'FontSize',12,'FontName','Times')
  set(gca,'YTick',0:2:10,'FontSize',12,'FontName','Times')
  %set(gca,'FontSize',20,'FontName','Times')
  xlim([4 15])
  ylim([0 10]) 
 
 
 subplot(2,3,6)
 plot(u1(6,:),zcc(:),'-black','LineWidth',1.0);
 xlabel('$\it{u}/u_{*,hi}$','FontSize',16,'FontName','Times',...
       'Interpreter','latex');
 ylabel('$\it{z/z_{h}}$','FontSize',16,'FontName','Times','Interpreter','latex');


  set(gca,'XTick',[5:5:15],'FontSize',12,'FontName','Times')
  set(gca,'YTick',0:2:10,'FontSize',12,'FontName','Times')
  %set(gca,'FontSize',20,'FontName','Times')
  xlim([4 15])
  ylim([0 10]) 

end
