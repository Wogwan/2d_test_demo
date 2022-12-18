function [index] = get_index_of_state(Q_estim,state)
    for i=1:length(Q_estim)
       
        if ((state.vel >= Q_estim(i).ini_vel) && (state.vel <=  Q_estim(i).end_vel) &&...
                (state.pos >= Q_estim(i).ini_pos) && (state.pos <= Q_estim(i).end_pos))
            %we found the state
            index = i;
            break;
        end
    end
end

