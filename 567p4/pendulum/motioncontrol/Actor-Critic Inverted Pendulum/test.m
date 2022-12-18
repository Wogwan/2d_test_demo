clc
clear all

%define hyper parameters
S_sqr_initial = 20000;
df = 0.98; %discount factor
a_actor = 0.01;
b_actor = 5;%5-10
a_critic = 0.01;
b_critic = 5;%5-10
actor_th_err = 0.01; %theshold to split
actor_th_n = 100; %number of samples to split

critic_th_err = 0.1; %theshold to split
critic_th_n = 100; %number of samples to split

% Initial state-action space
N = 3; %number of divisions

critic_function = initialize_Q_function(N);
actor_function = initialize_actor_function(N,S_sqr_initial); %each state has only one action? a number or a range

n_episodes = 4000;
n_iterations = 500;

total_accumulated_reward = zeros(1,n_episodes);
x_reward = 1;
 
for n=1:n_episodes
    %Set initial state pos = pi, vel = 0
    state.pos = pi;
    state.vel = 0;
    fprintf('Learning - Episode %d\n',n);
    for i=1:n_iterations
        [actor_action,actor_index] = select_random_action_actor(actor_function,state);
        
        %execute action
        [new_state, reward] = perform_action(state, actor_action);
        
        %critic evaluates the action executed
        new_action = select_action_policy(actor_function,new_state);
        
        [Q_policy ,i_policy] =  getQValue(critic_function,new_state,new_action);
        q = reward + df * Q_policy;
        
        %get the Q value from the critic in state s and actor_action
        [Q, critic_index] = getQValue(critic_function,state,actor_action);
        
        %check whether the action executed had a better result
        if q > Q 
            a_target = actor_action;
        elseif q < Q
            a_target = arg_max_a(critic_function,state);
        end
        
        %update action in the actor  
        lr_actor = 1/(a_actor*actor_function(actor_index).n + b_actor); % calculate learning rate
        if(actor_function(actor_index).action + lr_actor*(a_target-actor_function(actor_index).action)>=-5&&  actor_function(actor_index).action + lr_actor*(a_target-actor_function(actor_index).action)<=5)
            actor_function(actor_index).action =  actor_function(actor_index).action + lr_actor*(a_target-actor_function(actor_index).action); 
        end
        
        %update critic function value
        lr_critic = 1/(a_actor *critic_function(critic_index).n + b_actor); % calculate learning rate
        critic_function(critic_index).Q = critic_function(critic_index).Q + lr_critic * (q - critic_function(critic_index).Q);
       
        %update number of samples in the part
        actor_function(actor_index).n = actor_function(actor_index).n + 1;
        critic_function(critic_index).n = critic_function(critic_index).n + 1;
        
        %update variance
        actor_function(actor_index).S_sqr = actor_function(actor_index).S_sqr + lr_actor * (( a_target - actor_function(actor_index).action)^2 -  actor_function(actor_index).S_sqr );

               
        %update the state-action space by dividing it 
        actor_function = update_VR_actor(actor_function,actor_index,actor_th_err,actor_th_n,a_target,state);
        critic_function = update_VR_critic(critic_function,q,critic_index,critic_th_err,critic_th_n,state);

        %update current state
        state = new_state;
    end
     %test every 10 episodes
    if(mod(n,10)==0)
        fprintf('Testing - Episode %d\n',n);
        %%test using the policy
        state.pos = pi;
        state.vel = 0;
        for i=1:n_iterations
           %fprintf('Testing - Episode %d - Iteration %d\n ',n,  i);
           [action,index] = select_action_policy(actor_function,state);
           [new_state, reward] = perform_action(state, action);
           total_accumulated_reward(n) = total_accumulated_reward(n) + reward;
           %update current state
           state = new_state;
        end
    end
end


