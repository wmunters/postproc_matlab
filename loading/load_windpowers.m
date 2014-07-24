function [power_matrix] = load_windpowers
%loadwindpowers Loads power extracted from ABL by wind turbines
%   Function to load windpower from an existing file called Windpower.dat to a certain output matrix
%   The first two indices correspond to the i and j indices of the turbine array itself.
WP = load('Windpower.dat');
WP = -WP;
Nx = 8
Ny = 6

Nturb = Nx*Ny

power = WP(:,3:2:end);

for k=1:Nturb
    i = fix((k-1)/Ny)+1;
    j = k - (i-1)*Ny;
    
    power_matrix(i,j,:) = power(:,k);
    
    
end
