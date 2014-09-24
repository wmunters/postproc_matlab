function makepertmovie_zplane(Nx, Ny, Lx, Ly)
istart = 1.5;
istep = 0.01;
istop = 5.0;

turbines = false;
if(turbines)
    windfarm = load('windfarm.setup');
end

xgrid = makegrid(Lx,Nx);
ygrid = makegrid(Ly,Ny);
%close all
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
    
    
    clf
    hold on 
    filename = strcat('u_zplane_k008_t_',timename,'.dat');
    u = load(filename); u = reshape(u',[numel(u) 1]); u = reshape(u,[Nx Ny]);   
    
    %filename = strcat('u_yplane_j010_t_',timename,'.dat');
    %u = load(filename); u = reshape(u',[numel(u) 1]); u = reshape(u,[Nx Nz]);
%    subplot(2,1,1)
    title(strcat('U at t=',num2str(1000*i)))
    pcolor(xgrid, ygrid, u'); shading interp; axis equal; axis tight;
    caxis([13 20]); colorbar
    if(turbines)
        plot_turbines_topview(windfarm);
    end
    drawnow
   % pause
%    ylim([4 5]);xlim([1 2]);
%    
%    subplot(2,1,2)
%    pcolor(xgrid,ygrid,u'); shading interp; 
%    caxis([10 20]);
%    axis equal; axis tight; xlim([10 20]);
    
%     saveas(gca,strcat(filename,'.png'));
    
    
end
