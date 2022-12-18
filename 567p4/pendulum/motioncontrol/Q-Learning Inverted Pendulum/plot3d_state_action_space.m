function [] = plot3d_state_action_space(Q_func)
figure
hold on
for i=1:length(Q_func)
    plot3(Q_func(i).ini_vel,Q_func(i).end_vel,Q_func(i).ini_action,'-or')

axis([-pi-1 pi+1 -6 6 -2*pi-1 2*pi+1])
xlabel('Position') 
ylabel('Torque') 
zlabel('Velocity')
end

