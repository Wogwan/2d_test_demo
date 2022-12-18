function [Q_estim] = initialize_actor_function(N,S_sqr_initial)

max_vel = 2*pi;
min_vel = -2*pi;
max_pos = pi;
min_pos = -pi;
max_torque = 5;
min_torque = -5;

vel_increment = (abs(max_vel)+abs(min_vel))/N;
pos_increment = (abs(max_pos)+abs(min_pos))/N;
k=1;

for i=1:N
    for j=1:N
            Q_estim(k).ini_vel = min_vel + (i-1) * vel_increment;
            Q_estim(k).end_vel = min_vel + (i) * vel_increment;
            Q_estim(k).ini_pos =  min_pos + (j-1) * pos_increment;
            Q_estim(k).end_pos =  min_pos + (j) * pos_increment;
            Q_estim(k).action =  0;
            Q_estim(k).S_sqr = S_sqr_initial; %variance^2 - 22801
            Q_estim(k).n=1; 
            k=k+1;
    end
end
end

