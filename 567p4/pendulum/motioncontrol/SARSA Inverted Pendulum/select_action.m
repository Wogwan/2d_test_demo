function [action,index] = select_action(Q_estim, state, epsilon)
    
    %find all possible actions to the current state
    random_number = rand;
    if random_number>epsilon
        %exploit
        Qmax = -inf;
        selected_index = -1;

        for i=1:length(Q_estim)
            if ((state.vel >= Q_estim(i).ini_vel) && (state.vel <=  Q_estim(i).end_vel) &&...
                (state.pos >= Q_estim(i).ini_pos) && (state.pos <= Q_estim(i).end_pos))
                %we found the part for this state
                if(Q_estim(i).Q > Qmax)
                    selected_index = i;
                end
            end
        end
        action = ( Q_estim(selected_index).end_action + Q_estim(selected_index).ini_action )/2;
        index = selected_index;
    else
        %explore
        possible_actions = zeros(1,length(Q_estim));
        len=0;
        for i=1:length(Q_estim)
            if ((state.vel >= Q_estim(i).ini_vel) && (state.vel <=  Q_estim(i).end_vel) &&...
                    (state.pos >= Q_estim(i).ini_pos) && (state.pos <= Q_estim(i).end_pos))
                %we found the part for this state
                possible_actions(len+1) = i;
                len = len+1;
            end
        end
        index_random = possible_actions(randi(len));
        action =  Q_estim(index_random).ini_action + (Q_estim(index_random).end_action-Q_estim(index_random).ini_action)*rand();
        index = index_random;
    end        
    %action = ( Q_estim(index).end_action + Q_estim(index).ini_action )/2; %using mean approach to select the action
end

