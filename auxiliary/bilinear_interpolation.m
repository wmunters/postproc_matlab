function [u] = bilinear_interpolation(field, xmesh, ymesh, x, y)

	dx = xmesh(2) - xmesh(1);
	dy = ymesh(2) - ymesh(1);

	i_b = floor(x/dx) + 1;
	i_t = ceil(x/dx) + 1;
	j_b = floor(y/dy) + 1;
	j_t = ceil(y/dy) + 1;

    if(i_t > length(xmesh))
        i_t = i_t-length(xmesh);
    end
    if(j_t > length(ymesh))
        j_t = j_t -length(ymesh);        
    end
    
	x0 = xmesh(i_b);
	x1 = xmesh(i_t);
	y0 = ymesh(j_b);
	y1 = ymesh(j_t);

% watch out here, misschien u transposed doorgeven in plaats van u ...	
	uA = field(i_b,j_b);
	uB = field(i_b,j_t);
	uC = field(i_t,j_t);
	uD = field(i_t,j_b);
	
	if(x1 == x0) 
		cx = 0.0;
	else
		cx = (x-x0)/(x1-x0);
%         if(cx~=1 || cx~=0)
%             disp('Interpolating x')
%             cx
%         end
	end
	
	if(y1 == y0)
		cy = 0.0;
	else
		cy = (y-y0)/(y1-y0);
%         if(cy~=1 || cy~=0) 
%             disp('Interpolating y')
%             cy
%         end
	end
	
	% Interpolation in x direction 
	u_AD = uA + cx*(uD - uA);
	u_BC = uB + cx*(uC - uB);
	
	% Interpolation in y direction
	u = u_AD + cy*(u_BC - u_AD);

%     if(u>max([uA uB uC uD]) || u< min([uA uB uC uD]))
%         disp('INTERP')
%     end
end