function [k,uh_avg] = avg_spect(Nx, Ny, Lx, Ly, istart, istop, istep)

uh_avg = zeros(Nx,Ny/2-1);
counter = 0;

for i=istart:istep:istop
    i
    counter = counter+1;
    timename = adapt_timename(i,4);
    filename = strcat('u_zplane_k008_t_',timename,'.dat');

    u = load_plane(filename, Nx, Ny); 

    [k, uh] = pfft(u,Ly);

    uh_avg = uh_avg + uh;

end

uh_avg = uh_avg/counter;
