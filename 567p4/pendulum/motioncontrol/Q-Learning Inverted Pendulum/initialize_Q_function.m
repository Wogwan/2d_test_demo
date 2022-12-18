function [Q_estim] = initialize_Q_function(N,S_sqr_initial)

max_vel = 2*pi;
min_vel = -2*pi;
max_pos = pi;
min_pos = -pi;
max_torque = 5;
min_torque = -5;

torque_increment  = (abs(max_torque)+abs(min_torque))/N;
vel_increment = (abs(max_vel)+abs(min_vel))/N;
pos_increment = (abs(max_pos)+abs(min_pos))/N;
k=1;

for i=1:N
    for j=1:N
        for x=1:N 
            Q_estim(k).ini_vel = min_vel + (i-1) * vel_increment;
            Q_estim(k).end_vel = min_vel + (i) * vel_increment;
            Q_estim(k).ini_pos =  min_pos + (j-1) * pos_increment;
            Q_estim(k).end_pos =  min_pos + (j) * pos_increment;
            Q_estim(k).ini_action =  min_torque + (x-1) * torque_increment;
            Q_estim(k).end_action =  min_torque + (x) * torque_increment;

            Q_estim(k).n = 1; %should be 0 or 1??
            Q_estim(k).S_sqr = S_sqr_initial; %variance^2 - 22801
            Q_estim(k).Q = 0; %Q mean value
            k=k+1;
        end
    end
end
end

