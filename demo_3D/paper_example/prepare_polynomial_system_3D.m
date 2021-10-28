%%
clear all;
dbstop if error
format long

%%
syms x1 x2
dttr = 0.1;                                               % Recording step size for taining data (default = 0.3)
Ttr = 20;                                                 % Simulation time for training per starting point (default = 3)
noise = 0.01;                                             % Obervation noise
sn = noise*[1 1 1]';                                        % Observation noise (default = 1e-1)

%% Chebyshev interpolants value
tic
sz = 3;
deg = 4;
poly_deg = 4;
it = 400;

%% The first dimension
f1_p = -x1+x2;
f1_np = 0;
f1 = f1_p + f1_np;
y = chebfun(char(f1_np),[-sz,sz],'splitting','on'); % Modified here with different non-polynomial g
y_deg = minimax(y,deg); c_deg = chebcoeffs(y_deg);
T = chebyshevT([0:deg],x1);
f1_mid = vpa(T*c_deg);
x_change = x1/sz;
f1_appro_data = subs(f1_mid,x1,x_change);
f1_appro = f1_appro_data + f1_p;

%% The second dimension
f2_p = -x1-x2+x1*x2;
f2_np = x1*cos(x1);
f2 = f2_p + f2_np;
y = chebfun(char(f2_np),[-sz,sz],'splitting','on'); % Modified here with different non-polynomial g
y_deg = minimax(y,deg); c_deg = chebcoeffs(y_deg);
if length(c_deg)~=(deg+1)
    c_deg = [c_deg;0];
end
T = chebyshevT([0:deg],x1);
f2_mid = vpa(T*c_deg);
x_change = x1/sz;
f2_appro_data = subs(f2_mid,x1,x_change);
f2_appro = f2_appro_data + f2_p;
toc

%%
% Set parameters
dXtr_0 = [];                                              % Collecting training model for learning the difference in real and approximated dynamic systems
Xtr_0 = [];                                               % Collecting the Xtr_1 = [x1 x2] data by ode45 with setting
x0tr = [-0.5 -0.4; 0.2 -0.4];

E = 2;                                                    % Dimensions of the state space
dynt = @(t,x) dyn_controller_paper(0,x);                                 % dynamical system to be learned
dyn = @(x) dynt(0,x);                                     % time independent handle to dynamics

% Generate Training data
for i = 1:length(x0tr(1,:))
    if i == 1
        [t,xtr] = ode45(dynt,0:dttr:Ttr,x0tr(:,i)'); xtr = xtr';                                     % obtain the trajectories from the ODE45 with given sample time and given sample time-step
        x = xtr(:,1:end-1)';
        dtr = (xtr(:,2:end)-xtr(:,1:end-1))/dttr;
        %         noise_over_measurement = mvnrnd(zeros(E,1),diag(sn.^2),ntr)';                                                 % Obtain the xdot not directly, but with approximated differential method
        %         real_dtr = dtr + noise_over_measurement;
        %         real_dtr = dtr;
        d1_error = double(subs(f1_appro,{x1,x2},{x(:,1),x(:,2)}));
        d2_error = double(subs(f2_appro,{x1,x2},{x(:,1),x(:,2)}));
        y1 = dtr(1,:)'-d1_error;
        y2 = dtr(2,:)'-d2_error;
    elseif i == 2
        [t,xtr] = ode45(dynt,0:dttr:Ttr,x0tr(:,i)'); xtr = xtr';
        xtest = xtr(:,1:end-1)';
        dtr = (xtr(:,2:end)-xtr(:,1:end-1))/dttr;
        %         noise_over_measurement = mvnrnd(zeros(E,1),diag(sn.^2),ntr)';                                                 % Obtain the xdot not directly, but with approximated differential method
        %         real_dtr = dtr + noise_over_measurement;
        %         real_dtr = dtr;
        d1_error = double(subs(f1_appro,{x1,x2},{x(:,1),x(:,2)}));
        d2_error = double(subs(f2_appro,{x1,x2},{x(:,1),x(:,2)}));
        y1_test = dtr(1,:)'-d1_error;
        y2_test = dtr(2,:)'-d2_error;
    end
end

%%
dxdt1 = 0;
X = [1];
[mean1,hyp1,delta1,rmse1] = gpr_xdot1(x,y1,xtest,y1_test,it,noise,poly_deg);            % GP learning for the xdot2 dynamic systems
pvar x1 x2
X1 = p2s(monomials(x1,(0:poly_deg)));
X2 = p2s(monomials(x2,(0:poly_deg)));
for i = 2:length(X1)
    X = [X; X1(i); X2(i)];
end
for i = 1:length(X)
    dxdt1 = vpa(dxdt1 + X(i)*mean1(i));
end

%%
dxdt2 = 0;
X = [1];
[mean2,hyp2,delta2,rmse2] = gpr_xdot2(x,y2,xtest,y2_test,it,noise,poly_deg);            % GP learning for the xdot2 dynamic systems
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
figure(3);clf;hold on;
%%
subplot(211);hold on;
f1_learn = f1_appro + dxdt1;
syms x1 x2;
y1_learn = double(subs(f1_learn,{x1,x2},{xtest(:,1),xtest(:,2)}));
plot3(xtest(:,1),xtest(:,2),y1_learn(:,1),'b*'); hold on;
y1 = double(subs(f2,{x1,x2},{xtest(:,1),xtest(:,2)}));
plot3(xtest(:,1),xtest(:,2),y1(:,1),'ro'); hold on; view(30,40)
title(['rsme1 = ', num2str(rmse1)])

%%
subplot(212);hold on;
f2_learn = f2_appro + dxdt2;
syms x1 x2;
y2_learn = double(subs(f2_learn,{x1,x2},{xtest(:,1),xtest(:,2)}));
plot3(xtest(:,1),xtest(:,2),y2_learn(:,1),'b*'); hold on;
y2 = double(subs(f2,{x1,x2},{xtest(:,1),xtest(:,2)}));
plot3(xtest(:,1),xtest(:,2),y2(:,1),'ro'); hold on; view(30,40)
title(['rsme2 = ', num2str(rmse2)])
f_output = [f1_learn;f2_learn];

% g_chen = x1^2*x2-0.000006835767242-0.499957335901158*x1-0.125778690105319*x1^2+0.219390279834309*x1^3+0.220632156887092*x1^4+0.716967859965359*x1^5;
% g_chen = g_appro-0.000012375514062-0.500513008926802*x1+0.000542745343764*x2-0.119189589668064*x1^2-0.042173724594621*x1*x2+0.040714585039343*x2^2+0.305748528980397*x1^3+0.536554408293197*x1^2*x2+0.091448317009191*x1*x2^2-0.100083169093395*x2^3;
% g_chen = x1^2*x2-0.000007031474139-0.500009736044418*x1+0.126057255453549*x1^2+0.103008552092872*x1^3-0.168853172050168 *x1^4+0.264932312148793*x1^5;
% y3 = double(subs(g_chen,{x1,x2},{xtest(:,1),xtest(:,2)}));
% plot3(xtest(:,1),xtest(:,2),y3(:,1),'m^'); hold on;
