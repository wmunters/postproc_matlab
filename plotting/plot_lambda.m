function [] = plot_lambda(Nx,Ny,Lx,Ly,isovalue)

    close all
    savefigs = true;

    Nz = 80;
    Lz = 1;

    xmesh = makegrid(Lx,Nx);
    ymesh = makegrid(Ly,Ny);
    zmesh = makegrid(Lz,Nz);
    [X,Y,Z] = meshgrid(xmesh,ymesh,zmesh);

    istart = 0.025;
    istep  = 0.025;
    istop  = 1.000;


    figure
    for i=istart:istep:istop

        timename = adapt_timename(i,7);
        timenamevel = adapt_timename(i,4);

        filename_lambda = strcat('lambda2_t=',timename);
        filename_veloc  = strcat('velocity_field_t_',timenamevel,'.dat');

        lambda = load_field(filename_lambda,Nx,Ny,Nz);
        u      = load_field(filename_veloc,Nx,Ny,Nz,3); u = squeeze(u(:,:,:,1));


        figure(1)
        isosurface(X,Y,Z,lambda,isovalue,u);

    end
