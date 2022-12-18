function [Q_estim_updated] = update_state_space(Q_estim,q,index,th_err,th_n,state,epsilon_init)

max_vel = 2*pi;
min_vel = -2*pi;
max_pos = pi;
min_pos = -pi;
max_torque = 5;
min_torque = -5;


    if((q-Q_estim(index).Q)^2 >= th_err && Q_estim(index).n > th_n)
       last_index = length(Q_estim)+1;
       %split
       vel_dim = normalize(Q_estim(index).end_vel,max_vel,min_vel) - normalize(Q_estim(index).ini_vel,max_vel,min_vel);
       pos_dim = normalize(Q_estim(index).end_pos,max_pos,min_pos) - normalize(Q_estim(index).ini_pos,max_pos,min_pos);
       torque_dim = normalize(Q_estim(index).end_action,max_torque,min_torque) - normalize(Q_estim(index).ini_action,max_torque,min_torque);
       
       vel = Q_estim(index).end_vel - Q_estim(index).ini_vel;
       pos = Q_estim(index).end_pos - Q_estim(index).ini_pos;
       torque = Q_estim(index).end_action - Q_estim(index).ini_action;
       
       if(torque_dim >= vel_dim...
               && torque_dim  >= pos_dim)
       %split torque
       Q_estim(last_index) = Q_estim(index);
       Q_estim(index).end_action = Q_estim(index).end_action - torque/2;
       Q_estim(last_index).ini_action = Q_estim(index).end_action;
       
       elseif(( vel_dim >= torque_dim &&...
                vel_dim >= pos_dim) )
        %split vel

       Q_estim(last_index) = Q_estim(index);
       Q_estim(index).end_vel = Q_estim(index).end_vel - vel/2;
       Q_estim(last_index).ini_vel = Q_estim(index).end_vel;
       
       elseif((pos_dim >=  torque_dim  && ...
               pos_dim >= vel_dim))
       %split pos
       Q_estim(last_index) = Q_estim(index);
       Q_estim(index).end_pos = Q_estim(index).end_pos - pos/2;
       Q_estim(last_index).ini_pos = Q_estim(index).end_pos;
       end
       %what are the variance and Q values ini?
       %check where is q after partition
       %Q = where q is
       %
       if(state.pos <= Q_estim(index).end_pos && state.pos >= Q_estim(index).ini_pos && ...
           state.vel <= Q_estim(index).end_vel && state.vel >= Q_estim(index).ini_vel)
           Q_estim(index).Q = q;
           Q_estim(last_index).n = 1;
           Q_estim(last_index).epsilon = epsilon_init;
       elseif (state.pos <= Q_estim(last_index).end_pos && state.pos >= Q_estim(last_index).ini_pos && ...
           state.vel <= Q_estim(last_index).end_vel && state.vel >= Q_estim(last_index).ini_vel)
           Q_estim(last_index).Q = q;
           Q_estim(index).n = 1;
           Q_estim(index).epsilon = epsilon_init;

       end

    end
    Q_estim_updated = Q_estim;
    
end

