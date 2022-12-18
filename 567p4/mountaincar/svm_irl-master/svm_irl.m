function svm_irl()
grid_size=3;
m=mdp;
m.grid_size=grid_size;
m=m.value_iter();
m=m.get_optimal_policy;
start_state=[1,1]; 
% trajectory=m.generate_trajectory(start_state,length);
macrocell=feature_generation(1,m)
display('expectation of features:');
mu_e=expectation_features(start_state,macrocell,m)


%% algo starts here
m1=mdp;m1.grid_size=grid_size;
policy={m.random_policy()};
m1.policy=policy{1};
% traj=m1.generate_trajectory(start_state,length)
mu=expectation_features(start_state,macrocell,m1);
t=inf;
n=size(mu,1);m=size(mu,2);
epochs=0;
w_list=[];

while (t>exp(-9))
    epochs=epochs+1
    cvx_begin
        variable w(n,1)
        variable t
        maximise(t)
        subject to
            norm(w,2)<=1
            w'*(mu-mu_e*ones(1,m))<=-t.*ones(1,m)
    cvx_end
    w=normc(w);
    set_reward(w,macrocell,m1)
    m1.new_states=set_reward(w,macrocell,m1);
    m1=m1.value_iter();
    m1=m1.get_optimal_policy;
    % traj=m1.generate_trajectory(start_state,length)
    mu=[mu,expectation_features(start_state,macrocell,m1)];
    n=size(mu,1);m=size(mu,2);
    m1.grid_plotter(m1.grid_size,m1.states./max(m1.states(:)));
    w_list=([w_list,w])   
end
mu_e
mu
t
end
function states=set_reward(w,macrocell,m1)
states=zeros(size(m1.states));
for i=1:m1.n_states
    [x,y]=ind2sub(size(m1.states),i);
    states(x,y)=transpose(w)*feature_indicator(macrocell,[x,y]);
end
m1.states=states;
end
function sum=expectation_features(start_state,macrocell,mdp_obj)
sum=0;i=0;
state=start_state;
    while (norm(mdp_obj.gamma^(i).*feature_indicator(macrocell,state))> exp(-17))
        sum=sum+mdp_obj.gamma^(i).*feature_indicator(macrocell,state);
        state=mdp_obj.move(state,mdp_obj.actions(mdp_obj.policy(num2str(state)),:));
        i=i+1;
    end
end
function macrocell=feature_generation(grid_size,obj)
% feature generation
x=[1:obj.grid_size];y=x;
[x,y] = meshgrid(x,y);
x=x(:);y=y(:);
region_size=1;
macrocell=[];
i=1;
j=0;
while (j+region_size)<=size(x,1)
    macrocell=cat(3,macrocell,[x(j+1:j+region_size) y(j+1:j+region_size)]);
    i=i+1;
    j=j+region_size;
end
if j<size(x,1)
    z=zeros(region_size,2);
    z(1:size(x,1)-j+1,:)=[x(j:size(x,1)) y(j:size(x,1))];
    macrocell=cat(3,macrocell,z);
end
end
function feature= feature_indicator(macrocell,state)
for i=1:size(macrocell,3)
    for j=1:size(macrocell,1)
        if macrocell(j,:,i)==state
            feature=zeros(size(macrocell,3),1);
            feature(i)=1;
            return
        end
    end
end
end
