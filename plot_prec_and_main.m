function plot_prec_and_main(xmeshbig,startxp,cutoffxp,ymeshbig,startyp,cutoffyp,var_prec, xmesh,startx,cutoffx,ymesh,starty,cutoffy,var_main, end_stream, alpha)

    % Precursor
    pcolor(xmeshbig(startxp:cutoffxp),ymeshbig(startyp:cutoffyp),var_prec(startxp:cutoffxp,startyp:cutoffyp)')

    % Main
    h=pcolor(xmesh(startx:cutoffx),ymesh(starty:cutoffy),var_main(startx:cutoffx,starty:cutoffy)');


    caxis manual
    cmin = min( min(min(var_prec)), min(min(var_main)));
    cmax = max( max(max(var_prec)), max(max(var_main)));
    caxis([cmin cmax])
    colorbar
    shading interp
    axis equal
    axis tight
    % Rotate Main
    rotate(h,[0 0 1],-alpha/pi*180,[0 0 0]);
    % Reposition Main
    xd = get(h,'XData');
    yd = get(h,'Ydata');
    % Position top in origin of streamwise recycle plane
    yd = yd + end_stream(2);
    set(h,'YData',yd);
    % Position left in origin of spanwise recycle plane
    xd = xd + end_stream(1);
    set(h,'XData',xd);
