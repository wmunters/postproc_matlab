function [field_comp] = conjugate_wavenumber_spectrum(field)
%conjugate_wavenumber_spectrum - function for completing the spectral field representation
%   Function will fill in the part of the spectrum which is not stored in the SP-wind code because we rely on the Hermitian symmetry of the Fourier transform of a real array, i.e. the spectral field does not contain the negative values in k_x
%
% Syntax:  [field_comp] = conjugate_wavenumber_spectrum(field)
%
% Inputs:
%    field      - spectral field without negative k_x
%
% Outputs:
%    field_comp - spectral field completed with negative k_x
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: load_BLfield

% Author: Wim Munters, KU Leuven, Dept. of Mechanical Engineering
% Email : wim.munters@kuleuven.be
% Date  : 27 aug 2014


    Nxh = size(field,1);
    Nyh = size(field,2)/2+1;
    Ny  = size(field,2);
    Nz  = size(field,3);
    
    field_comp = zeros(Nxh-2,Ny,Nz);
    % ky=0
    field_comp(:,1,:) = conj(flipud(squeeze(field(2:end-1,1,:))));
    % ky>0
    for k=1:Nz
    field_comp(:,2:Nyh-1,k) = conj(fliplr(flipud(squeeze(field(2:end-1,Nyh+1:end,k)))));
    
    % ky<0
    field_comp(:,Nyh+1:end,k) = conj(fliplr(flipud(squeeze(field(2:end-1,2:Nyh-1,k)))));
    end
