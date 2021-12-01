%%
clear; dbstop if error; format long; tic;
%%
syms x1 x2 x3
dttr = 0.05;                                               % Recording step size for taining data (default = 0.3)
Ttr = 30;                                                 % Simulation time for training per starting point (default = 3)
noise = 1e-4;                                             % Obervation noise
sn = noise*[0 0 0]';                                        % Observation noise (default = 1e-1)
%% Chebyshev interpolants value
tic
sz = 5; deg = 4;
poly_deg = 4; it = 1000;
%% The first dimension
f1_p = -x1^2;
f1_np = -cos(x1)^2*sin(x1);
f1 = f1_p + f1_np;
y1 = chebfun(char(f1_np),[-sz,sz],'splitting','on'); % Modified here with different non-polynomial g
[y_deg_1, err1] = minimax(y1,deg); c_deg_1 = chebcoeffs(y_deg_1);
T = chebyshevT([0:deg],x1);
if length(c_deg_1)~=length(T)
    c_deg_1 = [c_deg_1;0];
end
f1_mid = vpa(T*c_deg_1);
x_change = x1/sz;
f1_appro_data = subs(f1_mid,x1,x_change);
f1_appro = f1_appro_data + f1_p;
%% The second dimension
f2_p = -x2-x1^3*x2;
f2_np = 0;
f2 = f2_p + f2_np;
%% The second dimension
f3_p = -x1^2*x3;
f3_np = 1-sqrt(sqrt((exp(x3)*cos(x3))^2));
f3 = f3_p + f3_np;
y3 = chebfun(char(f3_np),[-sz,sz],'splitting','on'); % Modified here with different non-polynomial g
[y_deg_3, err3] = minimax(y3,deg); c_deg_3 = chebcoeffs(y_deg_3);
T = chebyshevT([0:deg],x3);
if length(c_deg_3)~=length(T)
    c_deg_3 = [c_deg_3;0];
end
f3_mid = vpa(T*c_deg_3);
x_change = x3/sz;
f3_appro_data = subs(f3_mid,x1,x_change);
f3_appro = f3_appro_data + f3_p;
%% Chebyshev remainder item
% [Cheb_Pn1, Cheb_Pn3] = sos_cheb(deg,sz);                  % Provide the non-polynomial function in sos_cheb.m
[M1,x_num1] = cheb_max_3D(f1_np,sz);                          % Check the maximum value in the given region of the approximated function
[M3,x_num2] = cheb_max_3D(f3_np,sz);                         % Check the maximum value in the given region of the approximated function
hold on; plot(x_num1,M1,'ko'); plot(x_num2,M3,'ko');
% rho1 = cheb_rho1(deg,sz); % rho2 = cheb_rho2(deg,sz);
rho1 = 1.33;                                                % Obtain rho value of 4-th order from cheb_rho.m
rho2 = 1.89;                                                 % Obtain rho value of 4-th order from cheb_rho.m
d1 = 4*M1*rho1.^(-deg)/(rho1-1.0);                            % Upper bound
d3 = 4*M3*rho2.^(-deg)/(rho2-1.0);                            % Upper bound
%%
f_input = [f1_appro;f2;f3_appro]+[d1;0;d3];
f_origin = [f1,f2,f3];
%% Plot
figure(1);clf;hold on;
subplot(121);hold on;
a1 = plot(y1,'k');hold on;
set(a1,'linewidth',3,'linestyle','-.');
a2 = plot(y_deg_1,'m--');hold on;
set(a2,'linewidth',3,'linestyle',':');
set(a2,'color',[247 77 77]/255);
set(gca,'xtick',[-3,0,3]);
set(gca,'FontSize',20,'Fontname','Times');
set(gca,'Box','on');
ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 1.7;
xlabel('$x_1$','Interpreter','latex','Fontsize',22,'Fontname','Times');
ylabel('$y$','Interpreter','latex','Fontsize',22,'Fontname','Times');
subplot(122);hold on;
a1 = plot(y3,'k');hold on;
set(a1,'linewidth',3,'linestyle','-.');
a2 = plot(y_deg_3,'m--');hold on;
set(a2,'linewidth',3,'linestyle',':');
set(a2,'color',[247 77 77]/255);
set(gca,'xtick',[-3,0,3]);
set(gca,'FontSize',20,'Fontname','Times');
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.7;
xlabel('$x_3$','Interpreter','latex','Fontsize',22,'Fontname','Times');
ylabel('$y$','Interpreter','latex','Fontsize',22,'Fontname','Times');
refreshdata; drawnow;
%% Set parameters
dXtr_0 = [];                                              % Collecting training model for learning the difference in real and approximated dynamic systems
Xtr_0 = [];                                               % Collecting the Xtr_1 = [x1 x2] data by ode45 with setting
x0tr = [-0.1 -0.1 -0.2; -0.1 -0.2 -0.1; 0.1 0.1 0.2];
ntr = floor(Ttr/dttr);
E = 3;                                                    % Dimensions of the state space
dynt = @(t,x) dyn_controller_paper_3D(0,x);                                 % dynamical system to be learned
dyn = @(x) dynt(0,x);                                     % time independent handle to dynamics

% Generate Training data
for i = 1:length(x0tr(1,:))
    if i == 1
        [t,xtr] = ode45(dynt,0:dttr:Ttr,x0tr(:,i)'); xtr = xtr';                                     % obtain the trajectories from the ODE45 with given sample time and given sample time-step
        x_initial = xtr(:,1:end-1)';
        dtr_initial = (xtr(:,2:end)-xtr(:,1:end-1))/dttr;
        %%
        %         noise_over_measurement_1 = mvnrnd(zeros(E,1),diag(sn.^2),length(dtr_initial))';                                                 % Obtain the xdot not directly, but with approximated differential method
        %         real_dtr = dtr_initial + noise_over_measurement_1;
        real_dtr = dtr_initial;
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
        %         noise_over_measurement = mvnrnd(zeros(E,1),diag(sn.^2),length(dtr))';                                                 % Obtain the xdot not directly, but with approximated differential method
        %         real_dtr_test = dtr + noise_over_measurement;
        real_dtr_test = dtr;
        %%
        %         d1_error = double(subs(f1_appro,{x1,x2},{x(:,1),x(:,2)}));
        %         d2_error = double(subs(f2_appro,{x1,x2},{x(:,1),x(:,2)}));
        %         y1_test = dtr(1,:)'-d1_error;
        %         y2_test = dtr(2,:)'-d2_error;
    end
end
%%
dXtr_x_mid_1 = []; dXtr_xtest_mid_1 = [];
dXtr_x_mid_3 = []; dXtr_xtest_mid_3 = [];
for num = 1:length(x_initial(:,1))
    dXtr_x_mid_1 = [dXtr_x_mid_1; vpa(subs(f1,[x1,x2,x3],[x_initial(num,1),x_initial(num,2),x_initial(num,3)]))];
    dXtr_x_mid_3 = [dXtr_x_mid_3; vpa(subs(f3,[x1,x2,x3],[x_initial(num,1),x_initial(num,2),x_initial(num,3)]))];
end
for num = 1:length(xtest_initial(:,1))
    dXtr_xtest_mid_1 = [dXtr_xtest_mid_1; vpa(subs(f1,[x1,x2,x3],[xtest_initial(num,1),xtest_initial(num,2),x_initial(num,3)]))];
    dXtr_xtest_mid_3 = [dXtr_xtest_mid_3; vpa(subs(f3,[x1,x2,x3],[xtest_initial(num,1),xtest_initial(num,2),x_initial(num,3)]))];
end
%%
y1_mid = reshape(double(dXtr_x_mid_1),1,[])';
y1_test_mid = reshape(double(dXtr_xtest_mid_1),1,[])';
y1_trainset = reshape(double(real_dtr(1,:)),1,[])';
y1_testset = reshape(double(real_dtr_test(1,:)),1,[])';
y1 = y1_mid - y1_trainset;
y1_test = y1_test_mid - y1_testset;
%%
y3_mid = reshape(double(dXtr_x_mid_3),1,[])';
y3_test_mid = reshape(double(dXtr_xtest_mid_3),1,[])';
y3_trainset = reshape(double(real_dtr(3,:)),1,[])';
y3_testset = reshape(double(real_dtr_test(3,:)),1,[])';
y3 = y3_mid - y3_trainset;
y3_test = y3_test_mid - y3_testset;

%%
dxdt1 = 0;
X = [1];
[mean1,hyp2,delta2,rmse2] = gpr_xdot1_3D(x_initial,y1,xtest_initial,y1_test,it,noise,poly_deg);            % GP learning for the xdot2 dynamic systems
pvar x1 x2 x3
X1 = p2s(monomials(x1,(0:poly_deg)));
X2 = p2s(monomials(x2,(0:poly_deg)));
X3 = p2s(monomials(x3,(0:poly_deg)));
for i = 2:length(X1)
    X = [X;X1(i);X2(i);X3(i)];
end
for i = 1:length(X)
    dxdt1 = vpa(dxdt1 + X(i)*mean1(i));
end
%%
dxdt3 = 0;
X = [1];
[mean3,hyp3,delta3,rmse3] = gpr_xdot3_3D(x_initial,y3,xtest_initial,y3_test,it,noise,poly_deg);            % GP learning for the xdot2 dynamic systems
pvar x1 x2 x3
X1 = p2s(monomials(x1,(0:poly_deg)));
X2 = p2s(monomials(x2,(0:poly_deg)));
X3 = p2s(monomials(x3,(0:poly_deg)));
for i = 2:length(X1)
    X = [X;X1(i);X2(i);X3(i)];
end
for i = 1:length(X)
    dxdt3 = vpa(dxdt3 + X(i)*mean3(i));
end

f_output = [f1_appro+dxdt1;f2;f3_appro+dxdt3];
toc
%%
% %%
% figure(667);clf;
% set(gcf,'Position',[100 100 1000 350]);
% [ax1] = subplot(121);hold on;
% pos1 = set([ax1],'Position',[0.130 0.11 0.3347 0.8150]);
% y2_real = noise_over_measurement(2,:)'+double(subs(f2,{x1,x2,x3},{xtest_initial(:,1),xtest_initial(:,2),xtest_initial(:,3)}));
% scatter3(xtest_initial(:,1),xtest_initial(:,2),xtest_initial(:,3),40,y2_real,'filled')    % draw the scatter plot
% ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
% % cb = colorbar;                                     % create and label the colorbar
% % cb.Label.String = 'Values of 3D function';
% xlabel('$x_1$','Interpreter','latex','Fontsize',15);
% ylabel('$x_2$','Interpreter','latex','Fontsize',15);
% zlabel('$x_3$','Interpreter','latex','Fontsize',15);
% set(gca,'FontSize',15,'Fontname','Times');
% [ax2] = subplot(122);hold on;
% pos2 = set([ax2],'Position',[0.530 0.11 0.3347 0.8150]);
% f2_learn = f2_appro + dxdt2;
% y2_learn = double(subs(f2_learn,{x1,x2,x3},{xtest_initial(:,1),xtest_initial(:,2),xtest_initial(:,3)}));
% scatter3(xtest_initial(:,1),xtest_initial(:,2),xtest_initial(:,3),40,y2_learn,'filled')    % draw the scatter plot
% ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
% cb = colorbar;                                     % create and label the colorbar
% cb.Label.String = 'Values of 3D function';
% xlabel('$x_1$','Interpreter','latex','Fontsize',15);
% ylabel('$x_2$','Interpreter','latex','Fontsize',15);
% zlabel('$x_3$','Interpreter','latex','Fontsize',15);
% set(gca,'FontSize',15,'Fontname','Times');
%
% %%
% figure(668);clf;
% set(gcf,'Position',[100 100 1000 350]);
% [ax1] = subplot(121);hold on;
% pos1 = set([ax1],'Position',[0.130 0.11 0.3347 0.8150]);
% y3_real = noise_over_measurement(3,:)'+double(subs(f3,{x1,x2,x3},{xtest_initial(:,1),xtest_initial(:,2),xtest_initial(:,3)}));
% scatter3(xtest_initial(:,1),xtest_initial(:,2),xtest_initial(:,3),40,y3_real,'filled')    % draw the scatter plot
% ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
% % cb = colorbar;                                     % create and label the colorbar
% % cb.Label.String = 'Values of 3D function';
% xlabel('$x_1$','Interpreter','latex','Fontsize',15);
% ylabel('$x_2$','Interpreter','latex','Fontsize',15);
% zlabel('$x_3$','Interpreter','latex','Fontsize',15);
% set(gca,'FontSize',15,'Fontname','Times');
% [ax2] = subplot(122);hold on;
% pos2 = set([ax2],'Position',[0.530 0.11 0.3347 0.8150]);
% f3_learn = f3_appro + dxdt3;
% y3_learn = double(subs(f3_learn,{x1,x2,x3},{xtest_initial(:,1),xtest_initial(:,2),xtest_initial(:,3)}));
% scatter3(xtest_initial(:,1),xtest_initial(:,2),xtest_initial(:,3),40,y3_learn,'filled')    % draw the scatter plot
% view(-31,14);
% xlabel('$x_1$','Interpreter','latex','Fontsize',15);
% ylabel('$x_2$','Interpreter','latex','Fontsize',15);
% zlabel('$x_3$','Interpreter','latex','Fontsize',15);
% set(gca,'FontSize',15,'Fontname','Times');
% ax = gca;
% ax.XDir = 'reverse';
% cb = colorbar;                                     % create and label the colorbar
% cb.Label.String = 'Values of 3D function';
% pos3 = set(cd,'Position',[0.895 0.14 0.022 0.72]); hold on;
%% The first dimension
% f1_p = x2+x3^2;
% f1_np = 0;
% f1 = f1_p + f1_np;
% f1_appro = f1;
%% The second dimension
% f2_p = x3-x1^2;
% f2_np = 1*sin(x1);
% f2 = f2_p + f2_np;
% y2 = chebfun(char(f2_np),[-sz,sz],'splitting','on'); % Modified here with different non-polynomial g
% [y_deg_2, err2] = minimax(y2,deg); c_deg_2 = chebcoeffs(y_deg_2);
% T = chebyshevT([0:deg],x1);
% if length(c_deg_2)~=length(T)
%     c_deg_2 = [c_deg_2;0];
% end
% f2_mid = vpa(T*c_deg_2);
% x_change = x1/sz;
% f2_appro_data = subs(f2_mid,x1,x_change);
% f2_appro = f2_appro_data + f2_p;
% % rho2 = cheb_rho(deg,sz,char(f2_np));
% % rho2 = 2.1;
% % f2_appro = f2_appro_data + f2_p + rho2;
%% The second dimension
% f3_p = -x1^2*x3;
% f3_np = 2*sin(x3);
% f3 = f3_p + f3_np;
% y3 = chebfun(char(f3_np),[-sz,sz],'splitting','on'); % Modified here with different non-polynomial g
% [y_deg_3, err3] = minimax(y3,deg); c_deg_3 = chebcoeffs(y_deg_3);
% T = chebyshevT([0:deg],x3);
% if length(c_deg_3)~=length(T)
%     c_deg_3 = [c_deg_3;0];
% end
% f3_mid = vpa(T*c_deg_3);
% x_change = x3/sz;
% f3_appro_data = subs(f3_mid,x1,x_change);
% f3_appro = f3_appro_data + f3_p;
% % rho3 = cheb_rho(deg,sz,char(f3_np));
% % rho3 = 3.51;
% % f3_appro = f3_appro_data + f3_p + rho3;
