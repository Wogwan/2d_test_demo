function [action,index ] = select_action_policy(Q_estim,state)
    
    %find all possible actions to the current state
    Qmax = -inf;
    index = 1;
    for i=1:length(Q_estim)
        if ((state.vel >= Q_estim(i).ini_vel) && (state.vel <=  Q_estim(i).end_vel) &&...
                (state.pos >= Q_estim(i).ini_pos) && (state.pos <= Q_estim(i).end_pos))
            %we found the part for this state
            
            %select the highest Q
            if(Q_estim(i).Q > Qmax)
                Qmax = Q_estim(i).Q;
                index = i;
            end
        end
    end
    action = ( Q_estim(index).end_action + Q_estim(index).ini_action ) /2; %using mean approach to select the action
end

