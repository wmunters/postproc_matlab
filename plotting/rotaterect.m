function rotaterect(g,alpha,x0,y0,Lprime,alpha0)

% Rotate Precursor
        rotate(g,[0 0 1],alpha/pi*180,[0 0 0]);
        % Reposition Precursor
        
        pos = get(g,'Position');
        xp = pos(1);
        yp = pos(2);
%         xp = get(g,'XData');
%         yp = get(g,'YData');
        % Vertical alignment
        shift_y = y0 - Lprime*sin(alpha0 + alpha);
        yp = yp + shift_y;
%         set(g,'YData',yp);
        % Horizontal alignment
        shift_x = x0 - Lprime*cos(alpha0 + alpha);
        xp = xp + shift_x;
%         set(g,'XData',xp);
        pos = [xp yp pos(3) pos(4)];
        set(g,'Position',pos);