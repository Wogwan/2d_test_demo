function [Q] = calculate_Q(Q_estim, state,action)
 %find all possible actions to the current state
    for i=1:length(Q_estim)
        if( state.vel >= Q_estim(i).ini_vel && state.vel <=  Q_estim(i).end_vel &&...
                state.pos >= Q_estim(i).ini_pos && state.pos <= Q_estim(i).end_pos &&...
            action >= Q_estim(i).ini_action &&action <=  Q_estim(i).end_action   )
        %we found a part for this state
        Q = Q_estim(i).Q;
        break;
        end
    end
end

