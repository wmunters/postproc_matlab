function [BLfield,fieldr] = load_BLfield(name, Nx, Ny, Nz, Nl)
%load_BLfield - function for loading in spectral field representation
%
% Syntax:  [BLfield, fieldr] = load_BLfield(name, Nx, Ny, Nz, Nl)
%
% Inputs:
%    name       - name of the file, usually just 'BL_field.dat'
%    Nx         - size of the field in x direction (as in NS.setup)
%    Ny         - size of the field in y direction
%    Nz         - size of the field in z direction
%
% Outputs:
%    BLfield    - structure of complex arrays of size Nx/2 + 1, Ny, Nz
%    fieldr     - structure of real arrays of size Nx, Ny, Nz
%
% Other m-files required: conjugate_wavenumber_spectrum.m
% Subfunctions: none
% MAT-files required: none
%
% See also: load_BLfieldstat, load_field

% Author: Wim Munters, KU Leuven, Dept. of Mechanical Engineering
% Email : wim.munters@kuleuven.be
% Date  : 27 aug 2014

%------------- BEGIN CODE --------------

    Nxh = Nx/2 + 1;
    Nyh = Ny/2 + 1;

    field_dumm = load(name);
    field_dummc = field_dumm(:,1) + 1i*field_dumm(:,2);
   
    % First of all, load in the data

    % u     
    start  = 1;
    amount = Nxh*Ny*Nz;
    BLfield.u = reshape(field_dummc(start:start+amount-1), [Nxh Ny Nz]);
    
    % v 
    start = start+amount;
    amount = Nxh*Ny*Nz;
    BLfield.v = reshape(field_dummc(start:start+amount-1), [Nxh Ny Nz]);

    % w 
    start = start+amount;
    amount = Nxh*Ny*(Nz-1);
    BLfield.w = reshape(field_dummc(start:start+amount-1), [Nxh Ny Nz-1]);

    % theta
    if(Nl>3)
        start = start+amount;
        amount = Nxh*Ny*Nz;
        BLfield.theta = reshape(field_dummc(start:start+amount-1), [Nxh Ny Nz]);
    end

    % Now transform it to the real domain

    % First, fill up the missing wavenumber content (negative k_x), knowing that the complex field is Hermitian
    u_comp = conjugate_wavenumber_spectrum(BLfield.u);
    v_comp = conjugate_wavenumber_spectrum(BLfield.v);
    w_comp = conjugate_wavenumber_spectrum(BLfield.w);
    
    BLfield.u = [BLfield.u; u_comp];
    BLfield.v = [BLfield.v; v_comp];
    BLfield.w = [BLfield.w; w_comp];
    
    scaling = Nx*Ny;

    % Then,  perform the inverse Fourier transform and store in the fieldr structure
    fieldr.u = abs(ifft2(BLfield.u))*scaling;
    fieldr.v = abs(ifft2(BLfield.v))*scaling;
    fieldr.w = abs(ifft2(BLfield.w))*scaling;

%------------- END OF CODE --------------
