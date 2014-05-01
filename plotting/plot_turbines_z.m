
load windfarm.setup
hold on
Lx = 4.5*pi;
Nx = 512;
Ly = 2*pi;
Ny = 256;
Lz = 1;
Nz = 80;
% 
 windfarm(:,1) = windfarm(:,1);
 windfarm(:,2) = windfarm(:,2); 

for i=1:size(windfarm,1)
    
    plot([windfarm(i,1) windfarm(i,1)],[windfarm(i,2)-0.06 windfarm(i,2)+0.06],'k','LineWidth',3)
       
    
end