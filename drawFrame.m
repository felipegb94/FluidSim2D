function [] = drawFrame(Pos, numParticles, boxWidth, boxHeight)
%drawFrame Plot particles at a specific time step
    
    totalParticles = size(Pos, 2);
    hold on;
    axis([-0.2 boxWidth+0.1 -0.2 boxHeight+0.1]);
    
    for i = 1:numParticles
       x = Pos(1,i);
       y = Pos(2,i);
       scatter(x,y,'MarkerEdgeColor',[0 .5 .5],...
                   'MarkerFaceColor',[0 .7 .7],...
                   'LineWidth',1.5);
    end
    
    for i = numParticles+1:totalParticles
       x = Pos(1,i);
       y = Pos(2,i);
       scatter(x,y,'MarkerEdgeColor',[0 0 0],...
                   'MarkerFaceColor',[0 0 0],...
                   'LineWidth',1.5);
    end
    
    hold off;
    

end

