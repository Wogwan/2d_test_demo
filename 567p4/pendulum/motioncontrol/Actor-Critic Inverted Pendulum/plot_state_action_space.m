function [] = plot_state_action_space(Q_estim,torque)

n = 1;
figure
hold on;
for i=1:length(Q_estim)
    if(torque>=Q_estim(i).ini_action && torque<= Q_estim(i).end_action)
      points(1,1:2)  = [Q_estim(i).ini_vel Q_estim(i).ini_pos ];
      points(2,1:2)  = [Q_estim(i).ini_vel Q_estim(i).end_pos ];
      points(3,1:2)  = [Q_estim(i).end_vel Q_estim(i).end_pos ];
      points(4,1:2)  = [Q_estim(i).end_vel Q_estim(i).ini_pos ];
      points(5,1:2)  = [Q_estim(i).ini_vel Q_estim(i).ini_pos ];
%       points = [ 1 1; 2 2; 3 3; 4 4];
      plot(points(:,1),points(:,2),'Color',[0,0,1]);
      drawnow;
    end
  
    
end
ylabel('Position [-pi, pi]');
xlabel('Velocity [-2pi, 2pi]'); 

end
