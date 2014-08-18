function plot_1Dspect(k,in)

    close all

    in_avg = mean(abs(in),2);

    figure
    loglog(k,in_avg);
    xlabel('k')
   
    hold on
    plot(k,100*k.^(-5/3),'k')
