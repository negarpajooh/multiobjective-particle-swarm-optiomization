function PlotCosts(pop,rep)

    pop_costs=[pop.Cost];
    plot(abs(pop_costs(1,:)),abs(pop_costs(2,:)),'b.','linewidth',2);
    
    hold on;
    
    rep_costs=[rep.Cost];
    plot(abs(rep_costs(1,:)),abs(rep_costs(2,:)),'r*','linewidth',4);
    
    xlabel('1st Objective');
    ylabel('2nd Objective');
    
    grid minor;
    
    hold off;

end