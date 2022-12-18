function [actor_function_updated] = update_VR_actor(actor_function,actor_index,actor_th_err,actor_th_n,a_target,state)

max_vel = 2*pi;
min_vel = -2*pi;
max_pos = pi;
min_pos = -pi;

    if((a_target-actor_function(actor_index).action)^2 >= actor_th_err && actor_function(actor_index).n > actor_th_n)
       last_index = length(actor_function)+1;
       %split
       vel_dim = normalize(actor_function(actor_index).end_vel,max_vel,min_vel) - normalize(actor_function(actor_index).ini_vel,max_vel,min_vel);
       pos_dim = normalize(actor_function(actor_index).end_pos,max_pos,min_pos) - normalize(actor_function(actor_index).ini_pos,max_pos,min_pos);
       
       vel = actor_function(actor_index).end_vel - actor_function(actor_index).ini_vel;
       pos = actor_function(actor_index).end_pos - actor_function(actor_index).ini_pos;

       if(vel_dim > pos_dim )
           %split vel
           actor_function(last_index) = actor_function(actor_index);
           actor_function(actor_index).end_vel = actor_function(actor_index).end_vel - vel/2;
           actor_function(last_index).ini_vel = actor_function(actor_index).end_vel;
       elseif( pos_dim >= vel_dim)
           %split pos
           actor_function(last_index) = actor_function(actor_index);
           actor_function(actor_index).end_pos = actor_function(actor_index).end_pos - pos/2;
           actor_function(last_index).ini_pos = actor_function(actor_index).end_pos;
       end
     
       if(state.pos <= actor_function(actor_index).end_pos && state.pos >= actor_function(actor_index).ini_pos && ...
           state.vel <= actor_function(actor_index).end_vel && state.vel >= actor_function(actor_index).ini_vel)
           actor_function(actor_index).action = a_target;
           actor_function(last_index).n = 1;

       elseif (state.pos <= actor_function(last_index).end_pos && state.pos >= actor_function(last_index).ini_pos && ...
           state.vel <= actor_function(last_index).end_vel && state.vel >= actor_function(last_index).ini_vel)
           actor_function(last_index).action = a_target;
           actor_function(actor_index).n = 1;
       end

    end
    actor_function_updated = actor_function;
    
end

