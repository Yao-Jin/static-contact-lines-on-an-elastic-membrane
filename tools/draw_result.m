

close all;
    
    plot(x,y(2:end-1),'-k','LineWidth',1);
    hold on;
    axis equal;
    xlim([-1,1]);
    xlabel('$x$','Interpreter','latex');
    ylim([-0.3,1]);
    ylabel('$y$','Interpreter','latex');
    para = linspace(theta-pi/2,3/2*pi-theta,200)';
    plot(ones(200,1)*circenterx+Radius*cos(para),ones(200,1)*circentery+Radius*sin(para),...
         '-r','LineWidth',1);
%     plot([xcl1;xcl2],[ycl1;ycl2],'*b','MarkerSize',10);
    legend('membrane','droplet surface');
    hold off;
    
    axes('position',[0.18,0.52,0.27,0.27]);

    
    plot(x,y(2:end-1),'-k','LineWidth',1);
    hold on;
    axis equal;
    xlim([-1,1]);
    ylim([-0.3,1]);
    para = linspace(theta-pi/2,3/2*pi-theta,200)';
    plot(ones(200,1)*circenterx+Radius*cos(para),ones(200,1)*circentery+Radius*sin(para),...
         '-r','LineWidth',1);
%     plot([xcl1;xcl2],[ycl1;ycl2],'*b','MarkerSize',10);
%     legend('membrane','droplet surface','contact points');
    hold off;

