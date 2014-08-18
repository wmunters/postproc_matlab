function makepertmovie_yplane(Nx,Nz,Lx,Lz)

istart = 0
istop  = 5
istep = 0.05;
xgrid = makegrid(Lx,Nx);
zgrid = makegrid(Lz,Nz);
%Nx = 512;
%Ny = 128;
%Nz = 80;

figure
for i=istart:istep:istop
    
    i
    timename = num2str(i);
    if(length(timename)==1) 
        timename = strcat(timename,'.0000');
    else
        while (length(timename)<6)
            timename = strcat(timename,'0');
        end
    end
    
    figure(1)
    
    filename = strcat('u_yplane_j010_t_',timename,'.dat');
    u = load(filename); u = reshape(u',[numel(u) 1]); u = reshape(u,[Nx Nz]);
    pcolor(xgrid,zgrid,u'); shading interp; axis equal; axis tight;
    caxis([10 25]);
    
   % saveas(gca,strcat(timename,'.png'));
    
    
end
