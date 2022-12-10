function f_output = p4_gp(sys)
%%
% Generate the polynomial approximated system of inverted pendulum system

%%
dbstop if error
format long
rng(1)

%%
plott = sys.plott;
syms x1 x2
dttr = 0.1;                                               % Recording step size for taining data (default = 0.3)
Ttr = 20;                                                 % Simulation time for training per starting point (default = 3)
noise = 1e-3;                                             % Obervation noise

%% Chebyshev interpolants value
tic
poly_deg = 4;
it = 600;
g = 10;
l = 10;
dim = 2;

%% The first dimension
f1_p = x2;
f1_np = 0;
f1 = f1_p + f1_np;

%% The second dimension
f2_p = 0;
f2_np = -g/l*sin(x1);
f2 = f2_p + f2_np;

%% Set parameters
dXtr_0 = [];                                              % Collecting training model for learning the difference in real and approximated dynamic systems
Xtr_0 = [];                                               % Collecting the Xtr_1 = [x1 x2] data by ode45 with setting
% x0tr = [-0.05 -0.05];
ntr = floor(Ttr/dttr);
dynt = @(t,x) dyn_controller_paper_1d(0,x);                                 % dynamical system to be learned
dyn = @(x) dynt(0,x);                                     % time independent handle to dynamics

% Generate Training data
x0tr = rand(1,2);
[t,xtr] = ode45(dynt,0:dttr:Ttr,x0tr);                                  % obtain the trajectories from the ODE45 with given sample time and given sample time-step
y = -g./l.*sin(xtr(:,1))+noise.*randn(length(xtr(:,1)),1);

% Generate Testing data
x0tr = rand(1,2);
[t_t,xtr_t] = ode45(dynt,0:dttr:Ttr,x0tr);                                  % obtain the trajectories from the ODE45 with given sample time and given sample time-step
y_t = -g./l.*sin(xtr_t(:,1))+noise.*randn(length(xtr_t(:,1)),1);


%%
dxdt2 = 0;
X = [1];
[mean2,hyp2,delta2,rmse2] = p4_gpr_xdot2(xtr(:,1),y,xtr_t(:,1),y_t,it,noise,poly_deg,plott);            % GP learning for the xdot2 dynamic systems
pvar x1 x2
X1 = p2s(monomials(x1,(0:poly_deg)));
X2 = p2s(monomials(x2,(0:poly_deg)));
for i = 2:length(X1)
    X = [X; X1(i)];
end
for i = 1:length(X)
    dxdt2 = vpa(dxdt2 + X(i)*mean2(i));
end
f2_learn = dxdt2;

%%
if plott==1
    figure(3);clf;hold on;

    syms x1 x2;
    
    y2_learn = [];
    for num = 1:length(xtr_t)
        num_0 = double(vpa(subs(f2_learn,[x1],[xtr_t(num,1)])));
        y2_learn = [y2_learn; num_0];
    end
    y2_learn = reshape(double(y2_learn),1,[])';
    a1 = plot3(xtr_t(:,1),xtr_t(:,2),y2_learn,'b*'); hold on;
    
    y2_o = [];
    for num = 1:length(xtr_t)
        num_0 = vpa(subs(f2,[x1],[xtr_t(num,1)]));
        y2_o = [y2_o; num_0];
    end
    y2_o = reshape(double(y2_o),1,[])';
    a2 = plot3(xtr_t(:,1),xtr_t(:,2),y2_o,'ro'); hold on; view(30,40)
    
    % legend([a1,a2],{'Learned Value','Exact Dynamics'}, 'Interpreter','latex','location','northeast');
    legend([a2,a1],{'Exact Dynamics','Learned Value'}, 'Interpreter','latex','location','northeast');
    
    view(235, 25);hold on;
    title('');
    xlabel('x1','Interpreter','latex','Fontsize',21);
    ylabel('x2','Interpreter','latex','Fontsize',21);
    zlabel('x3','Interpreter','latex','Fontsize',21);
    set(gca,'Box','on');
    ax = gca;
    ax.LineWidth = 1.2;
    set(gca,'FontSize',21,'Fontname','Times');
end

%%
f_output = [f1;f2_learn];
% magnify
end