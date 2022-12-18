function [new_state, reward] = perform_action(state, action)

% action = torque
% run simulation for 0.1 seconds
% return reward and the new state
new_state = simulator(state,action); %will run for 0.1 sec
%plot_pendulum(new_state.pos);

%%reward function
reward = -abs(new_state.pos);

end

