function [] = plot_state_action_space(Q_estim)

n = 1;
figure
hold on;
for i=1:length(Q_estim)
      x  = (Q_estim(i).ini_vel+Q_estim(i).end_vel)/2;
      y  = (Q_estim(i).ini_pos+Q_estim(i).end_pos)/2;
 
%       points = [ 1 1; 2 2; 3 3; 4 4];
      plot(x,y,'.','Color',[0,0,1]);
%       drawnow;
  
    
end
ylabel('Position [-pi, pi]');
xlabel('Velocity [-2pi, 2pi]'); 
axis([-2*pi 2*pi -pi pi]);
end
