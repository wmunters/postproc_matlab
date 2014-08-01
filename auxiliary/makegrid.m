function [ grid ] = makegrid( L, N )
%MAKEGRID Generate one-dimensional grid (finite-difference fashion)
%   Generate 1D grid resembling the grids used in the sp-wind code. The
%   grid start from 0 and ends at L*(1-1/N)

    grid = linspace(0,L*(1-1/N),N);


