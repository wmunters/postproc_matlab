function [u] = load_plane(name,N1,N2)

%    u = load(name);

    fileID = fopen(name);
    u = fread(fileID,'real*8');
    u = reshape(u',[size(u,1)*size(u,2) 1]);
    if(N1*N2 ~= size(u,1))
        disp('File longer than expected..')
    end
    u = reshape(u(1:N1*N2),[N1 N2]);
%    figure
%    pcolor(u')
%    shading interp
%    axis equal
%    axis tight

