
function frame= PlotCosts(pop,F)
global W SatFactor CostFactor
    popCosts=[pop.Cost];
    popCosts(2,:)=(SatFactor-popCosts(2,:)*SatFactor)* sum(W);
    popCosts(1,:)=popCosts(1,:)*CostFactor;
    FCosts = [F.Cost];
    FCosts(2,:) = (SatFactor-FCosts(2,:)*SatFactor)* sum(W);  
    FCosts(1,:) = FCosts(1,:)*CostFactor;  
    plot(popCosts(1,:),popCosts(2,:),'b*','MarkerSize',8);
    hold on
    plot(FCosts(1,:),FCosts(2,:),'r*','MarkerSize',8);
    %plot(Costs(1,:),-1*Costs(2,:)*sum(W),'r*','MarkerSize',8);
    xlabel('Cost');
    ylabel('Satisfaction');
    title('Non-dominated Solutions (F_{1}) : Frontier');
    grid on;
    set(gcf,'color','w'); % set figure background to white
    drawnow;
    frame = getframe(1);
    hold off

end