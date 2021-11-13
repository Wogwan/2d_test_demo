%%
clear;
dbstop if error
format long

%%
syms x1 x2
dttr = 0.1;                                               % Recording step size for taining data (default = 0.3)
Ttr = 30;                                                 % Simulation time for training per starting point (default = 3)
noise = 1e-3;                                             % Obervation noise
sn = noise*[0 1]';                                        % Observation noise (default = 1e-1)

%% Chebyshev interpolants value
tic
% sz = 3;
sz = 2;
deg = 4;
poly_deg = 4;
it = 400;

%% The first dimension
% f1_p = x1^2*x2^2-2*x1-x2;
f1_p = -x1+x2;
f1_np = 0;
f1 = f1_p + f1_np;
f1_appro = f1;

%%
% y = chebfun(char(f1_np),[-sz,sz],'splitting','on'); % Modified here with different non-polynomial g
% y_deg = minimax(y,deg); c_deg = chebcoeffs(y_deg);
% T = chebyshevT([0:deg],x1);
% f1_mid = vpa(T*c_deg);
% x_change = x1/sz;
% f1_appro_data = subs(f1_mid,x1,x_change);
% f1_appro = f1_appro_data + f1_p;
%% The second dimension
f2_p = x1^2*x2;
f2_np = 1-sqrt(sqrt((exp(x1)*cos(x1))^2));
% f2_np = 1*sin(x1);
f2 = f2_p + f2_np;
y = chebfun(char(f2_np),[-sz,sz],'splitting','on'); % Modified here with different non-polynomial g
y_deg = minimax(y,deg); c_deg = chebcoeffs(y_deg);
T = chebyshevT([0:deg],x1);
% if length(c_deg) ~= length(T)
%     c_deg = [c_deg;0];
% end
f2_mid = vpa(T*c_deg);
x_change = x1/sz;
f2_appro_data = subs(f2_mid,x1,x_change);
f2_appro = f2_appro_data + f2_p;
rho = 1.9;                                                % Obtain rho value of 4-th order from cheb_rho.m
[M,x_num] = cheb_max(f2_np,sz);   
d = 4*M*rho.^(-deg)/(rho-1.0);                            % Upper bound
%%
% f_input = [f1;f2_appro];  
f_input = [f1;f2_appro]+[0;d];  
% the hyper-parameters in f2_appro is passing in dyn_controller_paper_1d.m

%% Set parameters
dXtr_0 = [];                                              % Collecting training model for learning the difference in real and approximated dynamic systems
Xtr_0 = [];                                               % Collecting the Xtr_1 = [x1 x2] data by ode45 with setting
x0tr = [-0.05 -0.05; 0.05 -0.05];
ntr = floor(Ttr/dttr);
E = 2;                                                    % Dimensions of the state space
dynt = @(t,x) dyn_controller_paper_1d(0,x);                                 % dynamical system to be learned
dyn = @(x) dynt(0,x);                                     % time independent handle to dynamics

% Generate Training data
for i = 1:length(x0tr(1,:))
    if i == 1
        [t,xtr] = ode45(dynt,0:dttr:Ttr,x0tr(:,i)'); xtr = xtr';                                     % obtain the trajectories from the ODE45 with given sample time and given sample time-step
        x_initial = xtr(:,1:end-1)';
        dtr_initial = (xtr(:,2:end)-xtr(:,1:end-1))/dttr;
        %%
%         noise_over_measurement = mvnrnd(zeros(E,1),diag(sn.^2),ntr)'; 
        noise_over_measurement = 0;
        real_dtr = dtr_initial + noise_over_measurement;
        %         real_dtr = dtr;
        %%
        %         d1_error = double(subs(f1_appro,{x1,x2},{x(:,1),x(:,2)}));
        %         d2_error = double(subs(f2_appro,{x1,x2},{x(:,1),x(:,2)}));
        %         y1 = dtr(1,:)'-d1_error;
        %         y2 = dtr(2,:)'-d2_error;
    elseif i == 2
        [t,xtr] = ode45(dynt,0:dttr:Ttr,x0tr(:,i)'); xtr = xtr';
        xtest_initial = xtr(:,1:end-1)';
        dtr = (xtr(:,2:end)-xtr(:,1:end-1))/dttr;
        %%
%         noise_over_measurement = mvnrnd(zeros(E,1),diag(sn.^2),ntr)';                                                 % Obtain the xdot not directly, but with approximated differential method
        noise_over_measurement = 0;
        real_dtr_test = dtr + noise_over_measurement;
        %         real_dtr = dtr;
        %%
        %         d1_error = double(subs(f1_appro,{x1,x2},{x(:,1),x(:,2)}));
        %         d2_error = double(subs(f2_appro,{x1,x2},{x(:,1),x(:,2)}));
        %         y1_test = dtr(1,:)'-d1_error;
        %         y2_test = dtr(2,:)'-d2_error;
    end
end

dXtr_x_mid = [];
dXtr_xtest_mid = [];
f_origin = [f1,f2];
for num = 1:length(x_initial(:,1))
    %         num_O = Odyn2D(Xtr_1(num,:)');                        % Collecting the original trajectories from the given dataset x=[x1,x2]
    dXtr_x_mid = [dXtr_x_mid; vpa(subs(f2,[x1,x2],[x_initial(num,1),x_initial(num,2)]))];
    dXtr_xtest_mid = [dXtr_xtest_mid; vpa(subs(f2,[x1,x2],[xtest_initial(num,1),xtest_initial(num,2)]))];
end
y2_mid = reshape(double(dXtr_x_mid),1,[])';
y2_test_mid = reshape(double(dXtr_xtest_mid),1,[])';
y2_trainset = reshape(double(real_dtr(2,:)),1,[])';
y2_testset = reshape(double(real_dtr_test(2,:)),1,[])';

%%
y2 = y2_mid - y2_trainset;
y2_test = y2_test_mid - y2_testset;

%%
% dxdt1 = 0;
% X = [1];
% [mean1,hyp1,delta1,rmse1] = gpr_xdot1(x,y1,xtest,y1_test,it,noise,poly_deg);            % GP learning for the xdot2 dynamic systems
% pvar x1 x2
% X1 = p2s(monomials(x1,(0:poly_deg)));
% X2 = p2s(monomials(x2,(0:poly_deg)));
% for i = 2:length(X1)
%     X = [X; X1(i); X2(i)];
% end
% for i = 1:length(X)
%     dxdt1 = vpa(dxdt1 + X(i)*mean1(i));
% end

%%
dxdt2 = 0;
X = [1];
[mean2,hyp2,delta2,rmse2] = gpr_xdot2(x_initial,y2,xtest_initial,y2_test,it,noise,poly_deg);            % GP learning for the xdot2 dynamic systems
pvar x1 x2
X1 = p2s(monomials(x1,(0:poly_deg)));
X2 = p2s(monomials(x2,(0:poly_deg)));
for i = 2:length(X1)
    X = [X; X1(i); X2(i)];
end
for i = 1:length(X)
    dxdt2 = vpa(dxdt2 + X(i)*mean2(i));
end
%%
% subplot(211);hold on;
% f1_learn = f1_appro + dxdt1;
% syms x1 x2;
% y1_learn = double(subs(f1_learn,{x1,x2},{xtest(:,1),xtest(:,2)}));
% plot3(xtest(:,1),xtest(:,2),y1_learn(:,1),'b*'); hold on;
% y1 = double(subs(f2,{x1,x2},{xtest(:,1),xtest(:,2)}));
% plot3(xtest(:,1),xtest(:,2),y1(:,1),'ro'); hold on; view(30,40)
% title(['rsme1 = ', num2str(rmse1)])
%%
figure(3);clf;hold on;
f2_learn = f2_appro + dxdt2;
syms x1 x2;
y2_learn = double(subs(f2_learn,{x1,x2},{xtest_initial(:,1),xtest_initial(:,2)}));
a1 = plot3(xtest_initial(:,1),xtest_initial(:,2),y2_learn(:,1),'b*'); hold on;
y2 = double(subs(f2,{x1,x2},{xtest_initial(:,1),xtest_initial(:,2)}));
% a2 = plot3(xtest_initial(:,1),xtest_initial(:,2),y2(:,1)+noise_over_measurement(2,:)','ro'); hold on; view(30,40)
a2 = plot3(xtest_initial(:,1),xtest_initial(:,2),y2(:,1),'ro'); hold on; view(30,40)

f_output = [f1;f2_learn];
legend([a1,a2],{'Learned Value','Exact Dynamics'}, 'Interpreter','latex','location','northeast');

view(235, 25);hold on;
title('');
xlabel('$x_1$','Interpreter','latex','Fontsize',18);
ylabel('$x_2$','Interpreter','latex','Fontsize',18);
zlabel('$x_3$','Interpreter','latex','Fontsize',18);
set(gca,'xtick',[-0.4,-0.2,0,0.2,0.4]);
set(gca,'ytick',[-0.4,-0.2,0,0.2,0.4]);
set(gca,'ztick',[0,0.2,0.4,0.6,0.8]);
set(gca,'Box','on');
ax = gca;
ax.LineWidth = 1.2;
xlim([-0.5 0.1]); ylim([-0.5 0.1]); zlim([-0.1 1]);hold on;
set(gca,'FontSize',18,'Fontname','Times');