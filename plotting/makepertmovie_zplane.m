function makepertmovie_zplane(Nx, Ny, Lx, Ly)
istart = 0.05
istop  = 5
istep = 0.05;

xgrid = makegrid(Lx,Nx);
ygrid = makegrid(Ly,Ny);

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
    filename = strcat('u_zplane_k008_t_',timename,'.dat');
    u = load(filename); u = reshape(u',[numel(u) 1]); u = reshape(u,[Nx Ny]);   
    
    %filename = strcat('u_yplane_j010_t_',timename,'.dat');
    %u = load(filename); u = reshape(u',[numel(u) 1]); u = reshape(u,[Nx Nz]);
    subplot(2,1,1)
    pcolor(xgrid, ygrid, u'); shading interp; axis equal; axis tight;
    caxis([10 20]);
    
    subplot(2,1,2)
    pcolor(xgrid,ygrid,u'); shading interp; 
    caxis([10 20]);
    axis equal; axis tight; xlim([10 20]);
    
    saveas(gca,strcat(filename,'.png'));
    
    
end
