clc
clear all

%define hyper parameters
S_sqr_initial = 500;
df = 0.98; %discount factor
a = 0.1;
b = 10;%5-10
th_err = 0.1; %theshold to split
th_n = 50; %number of samples to split

% Initial state-action space
N = 2; %number of divisions
Q_estim = initialize_Q_function(N,S_sqr_initial);

n_episodes = 3000;
n_iterations = 500;

total_accumulated_reward = zeros(1,n_episodes);


 
for n=1:n_episodes
    %Set initial state pos = pi, vel = 0
    state.pos = pi;
    state.vel = 0;
    fprintf('Learning - Episode %d\n',n);
    for i=1:n_iterations
        %fprintf('Learning - Episode %d - Iteration %d\n ',n,  i);
        [action,index] = select_action(Q_estim,state);
        [new_state, reward] = perform_action(state, action);
        Qmax = calculate_Qmax(Q_estim,new_state);
        q = reward + df*Qmax;
        %update number of samples in the part
        Q_estim(index).n = Q_estim(index).n + 1;
        %update Qi value
        alpha = 1/(a *Q_estim(index).n + b); % calculate learning rate
        Q_estim(index).Q = Q_estim(index).Q + alpha * (q - Q_estim(index).Q);
        %update variance
        Q_estim(index).S_sqr = Q_estim(index).S_sqr + alpha * (( q - Q_estim(index).Q)^2 -  Q_estim(index).S_sqr );
    
        %update the state-action space by dividing it 
        Q_estim = update_state_space(Q_estim,q,index,th_err,th_n,state);
        
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
           [action,index] = select_action_policy(Q_estim,state);
           [new_state, reward] = perform_action(state, action);
           total_accumulated_reward(n) = total_accumulated_reward(n) + reward;
           %update current state
           state = new_state;
        end
    end
end


