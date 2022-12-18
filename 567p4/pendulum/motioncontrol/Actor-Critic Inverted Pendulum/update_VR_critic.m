function [critic_function_updated] = update_VR_critic(critic_function,q,critic_index,critic_th_err,critic_th_n,state)

max_vel = 2*pi;
min_vel = -2*pi;
max_pos = pi;
min_pos = -pi;
max_torque = 5;
min_torque = -5;


if((q-critic_function(critic_index).Q)^2 >= critic_th_err && critic_function(critic_index).n > critic_th_n)
       last_index = length(critic_function)+1;
       %split
       vel_dim = normalize(critic_function(critic_index).end_vel,max_vel,min_vel) - normalize(critic_function(critic_index).ini_vel,max_vel,min_vel);
       pos_dim = normalize(critic_function(critic_index).end_pos,max_pos,min_pos) - normalize(critic_function(critic_index).ini_pos,max_pos,min_pos);
       torque_dim = normalize(critic_function(critic_index).end_action,max_torque,min_torque) - normalize(critic_function(critic_index).ini_action,max_torque,min_torque);
       
       vel = critic_function(critic_index).end_vel - critic_function(critic_index).ini_vel;
       pos = critic_function(critic_index).end_pos - critic_function(critic_index).ini_pos;
       torque = critic_function(critic_index).end_action - critic_function(critic_index).ini_action;
       
       if(torque_dim >= vel_dim...
               && torque_dim  >= pos_dim)
       %split torque
       critic_function(last_index) = critic_function(critic_index);
       critic_function(critic_index).end_action = critic_function(critic_index).end_action - torque/2;
       critic_function(last_index).ini_action = critic_function(critic_index).end_action;
       
       elseif(( vel_dim >= torque_dim &&...
                vel_dim >= pos_dim) )
        %split vel

       critic_function(last_index) = critic_function(critic_index);
       critic_function(critic_index).end_vel = critic_function(critic_index).end_vel - vel/2;
       critic_function(last_index).ini_vel = critic_function(critic_index).end_vel;
       
       elseif((pos_dim >=  torque_dim  && ...
               pos_dim >= vel_dim))
       %split pos
       critic_function(last_index) = critic_function(critic_index);
       critic_function(critic_index).end_pos = critic_function(critic_index).end_pos - pos/2;
       critic_function(last_index).ini_pos = critic_function(critic_index).end_pos;
       end
       %what are the variance and Q values ini?
       %check where is q after partition
       %Q = where q is
       %
       if(state.pos <= critic_function(critic_index).end_pos && state.pos >= critic_function(critic_index).ini_pos && ...
           state.vel <= critic_function(critic_index).end_vel && state.vel >= critic_function(critic_index).ini_vel)
           critic_function(critic_index).Q = q;
           critic_function(last_index).n = 1;

       elseif (state.pos <= critic_function(last_index).end_pos && state.pos >= critic_function(last_index).ini_pos && ...
           state.vel <= critic_function(last_index).end_vel && state.vel >= critic_function(last_index).ini_vel)
           critic_function(last_index).Q = q;
           critic_function(critic_index).n = 1;
       end

    end
    critic_function_updated = critic_function;
    
end

