function plot_recycleplanes(hinge_stream, end_stream, hinge_span, end_span, end_corner, cutoffx, cutoffy, Nx, Ny)

    fracx = cutoffx/Nx;
fracy = cutoffy/Ny;
downright_small = end_stream*(1-fracx) + end_corner*fracx;
topleft_small = hinge_stream*(fracy) + end_stream*(1-fracy);
topright_small = topleft_small + downright_small - end_stream;

    % Full main domain
    plot([hinge_stream(1) end_stream(1)],[hinge_stream(2) end_stream(2)],'k','LineWidth',0.3)
    plot([hinge_span(1) end_span(1)],[hinge_span(2) end_span(2)],'k','LineWidth',0.3)
    plot([end_corner(1) end_span(1)],[end_corner(2) end_span(2)],'k','LineWidth',0.3)
    plot([end_corner(1) end_stream(1)],[end_corner(2) end_stream(2)],'k','LineWidth',0.3)
    % Add boundaries of smaller domain
    if(cutoffx ~= Nx || cutoffy ~= Ny)
        plot([downright_small(1) end_stream(1)],[downright_small(2) end_stream(2)],':k','LineWidth',0.5)
        plot([topleft_small(1) end_stream(1)],[topleft_small(2) end_stream(2)],':k','LineWidth',0.5)
        plot([topleft_small(1) topright_small(1)],[topleft_small(2) topright_small(2)],':k','LineWidth',0.5)
        plot([downright_small(1) topright_small(1)],[downright_small(2) topright_small(2)],':k','LineWidth',0.5)
    end
    
