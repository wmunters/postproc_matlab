function [k,out] = pfft(in,L)
% PFFT Function to do a 'proper' fft

% 'Proper' fft - Takes out defunct mode and mean. Does fftshift and allow for easy further postprocessing

    N = size(in,2);
    inhat = fft(in,[],2);  % Note the dimension here!
    inhat = inhat';
%    inhat = fftshift(inhat);

%    figure

    out = inhat(2:N/2,:);
    %k = linspace(0,1,N/2);
    k = [0:(N/2-1)]*2*pi/L;
    k = k(2:end);
%    pcolor(log(abs(out)));
    out = out.*conj(out);
    out = out';



    
