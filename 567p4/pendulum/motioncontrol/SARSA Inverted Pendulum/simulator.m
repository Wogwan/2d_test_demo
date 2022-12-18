function [new_state] = simulator(state,torque)

%inputs:
    %action = torque
    %state = [current_position,current_velocity] in rad or rad/s
%outputs:
    %nexts = [next_position, next_velocity] in rad or rad/s
    %reward = reward obtained in the new position

%constants 
ff = 0.01; %friction factor
m=1; %mass
l=1; %length
g=9.8; %gravity
dt=0.001;

%dynamics
current_position = state.pos;
current_velocity = state.vel;
%torque = [-5,5];

for i=1:100
    %acceleration
    current_aceleration = (-ff*current_velocity+m*g*l*sin(current_position) + torque)/(m*l^2);

    %velocity
    current_velocity = current_velocity + current_aceleration*dt;

    %position
    current_position = current_position + current_velocity*dt;

    if(current_position>=pi)
        current_position = current_position-2*pi;
    elseif(current_position<=-pi)
        current_position = current_position+2*pi;
    end

    %limit velocity
    if(current_velocity>2*pi)
        current_velocity=2*pi;
    elseif(current_velocity<-2*pi)
        current_velocity=-2*pi;
    end

end
%set new state
new_state.pos = current_position;
new_state.vel = current_velocity;
% disp(new_state)
% torque
end
