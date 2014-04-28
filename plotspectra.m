close all
clear all

sp  = load('BL_spec_stream_inst.dat');
spp = load('BL_spec_stream_inst_precursor.dat');

Nz = 80;
N = 128;
L = 2*pi;
k=[0:1:N-1];
k2 = k/L;
    [prefix,~] = strtok(path,';');
    prefix = strcat(prefix,'\');
    prefix=''


sp_uu = sp(1:Nz,2:end);
sp_vv = sp(Nz+1:2*Nz,2:end);
sp_ww = sp(2*Nz+1:3*Nz - 1,2:end);

spp_uu = spp(1:Nz,2:end);
spp_vv = spp(Nz+1:2*Nz,2:end);
spp_ww = spp(2*Nz+1:3*Nz - 1,2:end);

% figure

loglog(k2,sp_uu(8,:),'r')
hold on
loglog(k2,spp_uu(8,:),'b')
legend('Main','Precursor')
xlabel('k_x/L_x')
title('uu')
saveas(gcf,strcat(prefix,'uuspec.png'))
hold off


loglog(k2,sp_vv(8,:),'r')
hold on
loglog(k2,spp_vv(8,:),'b')
legend('Main','Precursor')
xlabel('k_x/L_x')
title('vv')
saveas(gcf,strcat(prefix,'vvspec.png'))
hold off



loglog(k2,sp_ww(8,:),'r')
hold on
loglog(k2,spp_ww(8,:),'b')
legend('Main','Precursor')
xlabel('k_x/L_x')
title('ww')
saveas(gcf,strcat(prefix,'wwspec.png'))
hold off