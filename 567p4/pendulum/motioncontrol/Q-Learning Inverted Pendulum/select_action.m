function [action,index ] = select_action(Q_estim,state)
    
    %find all possible actions to the current state
    Q_rand.Qmax = -inf;
    Q_rand.index = 1;
    for i=1:length(Q_estim)
       
        if ((state.vel >= Q_estim(i).ini_vel) && (state.vel <=  Q_estim(i).end_vel) &&...
                (state.pos >= Q_estim(i).ini_pos) && (state.pos <= Q_estim(i).end_pos))
            %we found the part for this state
            S = sqrt(Q_estim(i).S_sqr);
            Q = Q_estim(i).Q;
            n = Q_estim(i).n;
            Q_rand_calculated = normrnd(Q,S);
            if(Q_rand_calculated > Q_rand.Qmax)
                Q_rand.Qmax = Q_rand_calculated;
                Q_rand.index = i;
            end
        end
    end
    Qmax = Q_rand.Qmax;
    index = Q_rand.index;
%     action = ( Q_estim(index).end_action + Q_estim(index).ini_action )/2; %using mean approach to select the action
    action =  Q_estim(index).ini_action + (Q_estim(index).end_action-Q_estim(index).ini_action).*rand();
end

