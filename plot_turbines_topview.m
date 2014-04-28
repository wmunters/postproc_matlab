function plot_turbines_topview(windfarm)

 for i=1:size(windfarm,1)
    
     plot([windfarm(i,1) windfarm(i,1)],[windfarm(i,2)-windfarm(i,4) windfarm(i,2)+windfarm(i,4)],'k','LineWidth',2)
          
 end




end
