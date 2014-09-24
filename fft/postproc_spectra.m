function [uh_avg,k] = postproc_spectra(xlocs,uh_in, k_in)

    close all

    % Some HARDCODED stuff here...
    Nx = 2048;
    Ny = 128;
    Lx = 16*pi;
    Ly = pi;
    istart = 3;
    istep  = 0.01;
    istop  = 5;
    symbolslist = {'o'; 's'; '<'; '>'};
    colorslist = {'k'; 'b'; 'r'};
    lineslist = {'-'; '.-'};
    frmt = {'-ok'; '-<k'; '--sk'; '-ob'; '-<b'; '--sb'; '-or'; '-<r'; '--sr'; '-or'; '-<r'};
%    frmt = ['b'; 'r'; 'k'; 'g'; 'c'; 'm'];
%    frmt= {'-ok';'-sk';'-->k';':<k';':k';'-ob';'-sb';'--db';':xb';':b'};

    % xlocs are given in km of length
    xmesh = makegrid(Lx,Nx);
    for i=1:numel(xlocs)
        x_indices(i) = find(xmesh>=xlocs(i),1);
        end

    % First get the averaged spectra at every x location
    if(nargin==3) 
        uh_avg = uh_in;
        k = k_in;
        else
        [k, uh_avg] = avg_spect(Nx,Ny,Lx,Ly,istart,istop,istep);
        end

    % Then make a plot for the x locations given as input
    figure
    for i=1:numel(xlocs)
        
        %loglog(k,uh_avg(x_indices(i),:),frmt{i},'LineWidth',2,'MarkerFaceColor',frmt{i}(end)) 
        loglog(k,uh_avg(x_indices(i),:),frmt{i},'LineWidth',2)
        hold on

    end

    legendstrings = strcat(num2str(xlocs), ' km');

    legend(legendstrings)

    saveas(gca,'spectra.png');