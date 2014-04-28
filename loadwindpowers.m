function [power_matrix] = loadwindpowers
WP = load('Windpower.dat');
WP = -WP;
Nx = 6
Ny = 6

Nturb = Nx*Ny

power = WP(:,3:2:end);

for k=1:Nturb
    i = fix((k-1)/Ny)+1;
    j = k - (i-1)*Ny;
    
    power_matrix(i,j,:) = power(:,k);
    
    
end
