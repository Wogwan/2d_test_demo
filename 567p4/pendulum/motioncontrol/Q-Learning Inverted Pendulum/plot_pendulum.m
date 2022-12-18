function plot_pendulum(theta)
    x_pos = -cos(theta+pi/2);
    y_pos = sin(theta+pi/2);
    X1 = 0;
    Y1 = 0;
    r=0.7;
    X2 = r*x_pos;
    Y2 = r*y_pos;
    plot(X2,Y2,'o','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b') ;
    hold on;
    plot([X1; X2], [Y1; Y2],'b') ;
    axis([-1 1 -1 1]);
    drawnow;
    hold off;

end
