function [action,index ] = select_action_policy(actor_function,state)
    
    %find all possible actions to the current state
    for i=1:length(actor_function)
        if ((state.vel >= actor_function(i).ini_vel) && (state.vel <=  actor_function(i).end_vel) &&...
                (state.pos >= actor_function(i).ini_pos) && (state.pos <= actor_function(i).end_pos))
            %what happens if the pos or vel fall in a split?
            action = actor_function(i).action;
            index=i;
            break;
        end
    end
end

