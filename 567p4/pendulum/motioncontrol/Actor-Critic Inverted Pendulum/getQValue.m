function [Q,index]= getQValue(critic_function,state,action)
      
    for i=1:length(critic_function)
        if ((state.vel >= critic_function(i).ini_vel) && (state.vel <=  critic_function(i).end_vel) &&...
                (state.pos >= critic_function(i).ini_pos) && (state.pos <= critic_function(i).end_pos)&&...
                 (action >= critic_function(i).ini_action) && (action <= critic_function(i).end_action))
            %what happens if it's in the edge of the part??
            Q = critic_function(i).Q;
            index=i;
            break;
        end
    end
end

