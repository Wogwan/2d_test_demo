classdef mdp
   properties
        grid_size; 
        actions=[0,-1;1,0;0,1;-1,0;0,0];
        %action 1: go left,   2: go down(in matrix) go up in plot
        %action 3: go right,  4: go up( in matrix) go down in plot
        gamma=0.9;
        epsilon=0.01;
        V;
        policy=containers.Map;
        policy_vector;
        new_states=[];
        ploton=1;
   end
    properties(Dependent, SetAccess = public)
        reward;
        states;
        converge_factor;
        n_states;
    end
   methods
       function obj = set.states(obj,val)
          obj.new_states = val;
       end
        function val=get.converge_factor(obj)
            val=obj.epsilon*(1-obj.gamma)/obj.gamma;
        end
        function val=get.n_states(obj)
            val=size(obj.states(:),1);
        end
        function val=get.reward(obj)
            val=zeros(obj.grid_size+1,obj.grid_size+1);
            val(obj.grid_size,obj.grid_size)=1;
        end
        function val=get.states(obj)
           if isempty(obj.new_states)
               val=obj.reward(1:obj.grid_size,1:obj.grid_size);
           else
               val=obj.new_states;
           end
        end
        function obj=value_iter(obj)
            V_new=zeros(size(obj.states));
            obj.V=V_new;
            diff=(inf);
            epochs=0;
            while any(obj.converge_factor<diff(:))
%                 epochs=epochs+1
                for i=1:obj.n_states
                    %get indices of the state
                    [x,y]=ind2sub(size(obj.states),i);
                    max_val=(-inf);
                    for j=1:size(obj.actions,1)
                        ind=obj.move([x,y],obj.actions(j,:));
                        x_new=ind(1);y_new=ind(2);
                        if (max_val<obj.gamma*obj.V(x_new,y_new)+obj.states(x,y))
                            max_val=obj.gamma*obj.V(x_new,y_new)+obj.states(x,y);
                            best_action=j;
                        end
                    end
                    V_new(x,y)=max_val;
                end
                diff=V_new-obj.V;
                obj.V=V_new;
%                 display('value:')
%                 display(obj.V);
            end
            if obj.ploton
                plot_value_function(obj,obj.V,obj.reward);
                obj.get_optimal_policy()
            end
       end
       function out=move(obj,state,action)
           if (obj.get_legal(state+action,size(obj.states)))
               out=(state+action);
               return
           else
               out=state;
               return 
           end
       end
       function out=get_legal(obj,s_ind,state_size)
          if(s_ind<=state_size(1))&(0<s_ind)
              out=true; 
              return
          else
              out=false;
              return 
          end
       end
       function plot_value_function(obj,V,reward)
                figure(1);
                subplot(1,2,1)
                title('Value function update');
                obj.grid_plotter(obj.grid_size,V./max(V(:)));
                
                hold on;
                subplot(1,2,2)
                title('reward function');
                obj.grid_plotter(obj.grid_size,obj.states);
       end
       function obj=get_optimal_policy(obj,s)
            obj.policy_vector=zeros(size(obj.states));
            for i=1:obj.n_states
                %get indices of the state
                [x,y]=ind2sub(size(obj.states),i);
                max_val=(-inf);
                for j=1:size(obj.actions,1)
                    ind=obj.move([x,y],obj.actions(j,:));
                    x_new=ind(1);y_new=ind(2);
                    if (max_val<obj.V(x_new,y_new))
                        max_val=obj.V(x_new,y_new);
                        best_action=j;
                    end
                end
                obj.policy(num2str([x,y]))=best_action;
                obj.policy_vector(x,y)=best_action;
            end
            if obj.ploton
                figure(1); %%this plot has y of matrix becomes x and x of matrix becomes y
                subplot(1,2,1);
                hold on
                x=[1:obj.grid_size]+0.75;y=x;
                [x,y] = meshgrid(x,y);
                v=reshape(obj.actions(obj.policy_vector(:),1),size(x)); u=reshape(obj.actions(obj.policy_vector(:),2),size(x));
                quiver(x,y,u,v);
            end
       end
       function p1=random_policy(obj)
            p1=containers.Map;
            for i=1:obj.n_states
                %get indices of the state
                [x,y]=ind2sub(size(obj.states),i);
                rand_action=randi([1 size(obj.actions,1)],1,1);
                ind=obj.move([x,y],obj.actions(rand_action,:));
                x_new=ind(1);y_new=ind(2);
                p1(num2str([x,y]))=rand_action;
            end
       end
%        function traj=random_trajectory(obj,start_state, length)
%             traj=[];
%             state=start_state;
%             for i=1:length
%                action=randi([1 size(obj.actions,1)],1,1);
%                traj=[traj; state, action];
%                state=obj.move(state,obj.actions(action,:));
%             end
%         end
       function traj=generate_trajectory(obj,start_state,length)
           traj=[];
           state=start_state;
           for i=1:length
               obj.policy(num2str(state));
               traj=[traj; state, obj.policy(num2str(state))];
               state=obj.move(state,obj.actions(obj.policy(num2str(state)),:));
           end
       end
       function grid_plotter(obj,no_cells, matrix)
            colormap gray;
            cmap = (colormap); %get current colormap
        %     cmap=cmap([1 60],:); % set your range here
            colormap(cmap); % apply new colormap
            x = [1:no_cells+1];
            y = [1:no_cells+1];
            [X,Y] = meshgrid(x,y);
            Z=zeros(no_cells+1,no_cells+1);
            Z(1:no_cells,1:no_cells)=matrix;
            pcolor(X,Y,Z);
            drawnow;
       end
   end      
end