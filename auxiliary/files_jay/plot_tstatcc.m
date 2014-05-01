%clear all
%clc
A=importdata('BL_tstatcc.dat',',');
z=A(:,1); % read all the rows,1st colummn
Um=A(:,2); % read all the rows,2nd column
Vm=A(:,3);
Wm=A(:,4);
uu=A(:,5);
vv=A(:,6);
ww=A(:,7);
uv=A(:,8);
uw=A(:,9);
vw=A(:,10);

z_h = 0.1;  % hub height

%For multiplying two terms
%for m=1:80
 %aa(m,1)=z(m,1)*z(m,1);
 %aa(m,1)
 %end   
figure
semilogx(z/z_h,Um,'-black','LineWidth',2);

xlabel('$\it{z/z_{h}}$','FontSize',24,'FontName','Times','Interpreter','latex');
ylabel('$\it{\left\langle\overline{u}\right\rangle/u_{*,hi}}$','FontSize',20,'FontName','Times',...
       'Interpreter','latex')

 set(gca,'XTick',[0.1 1 10],'FontSize',18,'FontName','Times')
 set(gca,'YTick',4:2:16,'FontSize',18,'FontName','Times')
 set(gca,'FontSize',20,'FontName','Times')
 xlim([0.05 20])
 ylim([4 16])

line([1-1/2 1-1/2],[4 16],'LineStyle',':','Color','black','LineWidth',1.0)
line([1+1/2 1+1/2],[4 16],'LineStyle',':','Color','black','LineWidth',1.0)
   

figure
plot(-uw,z/z_h,'--black','LineWidth',1.5);
xlabel('$\it{\tau_{xz}/u_{*}^{2}}$','FontSize',20,'FontName','Times',...
       'Interpreter','latex')
ylabel('$\it{z/z_{h}}$','FontSize',20,'FontName','Times','Interpreter','latex');

 set(gca,'XTick',0.0:0.2:1.0,'FontSize',18,'FontName','Times')
 set(gca,'YTick',0:2:10,'FontSize',18,'FontName','Times')
 xlim([0.0 1.0])
 ylim([0 10])

line([0 1],[1-1/2 1-1/2],'LineStyle',':','Color','black','LineWidth',1.0)
line([0 1],[1+1/2 1+1/2],'LineStyle',':','Color','black','LineWidth',1.0)


