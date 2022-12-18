function [Qmax] = calculate_Qmax(Q_estim, state)
 %find all possible actions to the current state
    k=1;
    for i=1:length(Q_estim)
        if( state.vel >= Q_estim(i).ini_vel && state.vel <=  Q_estim(i).end_vel &&...
                state.pos >= Q_estim(i).ini_pos && state.pos <= Q_estim(i).end_pos)
        %we found a part for this state
        Q = Q_estim(i).Q;
        all_Q(k) = Q;
        k = k+1;
        end
    end
    Qmax = max(all_Q);
end

