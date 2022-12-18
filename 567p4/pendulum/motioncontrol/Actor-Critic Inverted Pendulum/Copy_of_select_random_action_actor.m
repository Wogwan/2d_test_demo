function [action,index] = select_random_action_actor(actor_function,state)
  
    %find all possible actions to the current state
    index = -1;
    possible_states = zeros(1,length(actor_function));
    possible_states_count=0;
    for i=1:length(actor_function)
       
        if ((state.vel >= actor_function(i).ini_vel) && (state.vel <=  actor_function(i).end_vel) &&...
                (state.pos >= actor_function(i).ini_pos) && (state.pos <= actor_function(i).end_pos))
            possible_states_count = possible_states_count+1;
            possible_states(possible_states_count) = i;
        end
    end
    %select a random state between the different possible states, this
    %happens when it falls in the edges of the parts
    %we found the part for this state 
    %select an action with some exploration
    index_rand = randi(possible_states_count);
    state_selected = actor_function(possible_states(index_rand));
    
    S = sqrt(state_selected.S_sqr);
    action = state_selected.action;
    a_rand_calculated = normrnd(action,S);
    action = a_rand_calculated;
    if(action>5)
        action = 5;
    elseif action<-5
        action=-5;
    end
    index = i;

end

