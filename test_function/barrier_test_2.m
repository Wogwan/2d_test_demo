% function [Bs, S] = barrier_test(k,f,gam_0,c)

pvar x1 x2 gamsym; 
%%
k = 4;
dom = 4;
domain = [-dom dom -dom dom];
gam_0 = 0.2;

c1 = gam_0 -(x1^2+x2^2); % Safe region qualified polynomial X_0
c2 = x1^2+x2^2-4;    % Unsafe region : X_u
c3 = x1^2+x2^2;       % For all [x x] : X
c4 = 4-c3;              % X \ X_u
c = [c1;c2;c3;c4];

f = [x2; 
    (1-x1^2)*x2-x1];

%%
% Define polynomial variables
      
x = [x1;x2];

[B,Bc] = polydecvar('b',monomials(x,0:k));
% [B,Bc] = sosdecvar('b',monomials(x,0:k));

[S1,S1_Q] = sosdecvar('S1_w',monomials(x,0:k)); % S2 sos decision variables
[S2,S2_Q] = sosdecvar('S2_w',monomials(x,0:k)); % S3 sos decision variables
[S3,S3_Q] = sosdecvar('S3_w',monomials(x,0:k)); % S4 sos decision variables
[S4,S4_Q] = sosdecvar('S4_w',monomials(x,0:k)); % S5 sos decision variables
[S5,S5_Q] = sosdecvar('S5_w',monomials(x,0:k)); % S6 sos decision variables

% Define SOS constraints
Mid = 0.5'*0.5*sum(diag(jacobian(jacobian(B,x),x)))*0.5;

pconstr(1) = S1>=0;
pconstr(2) = S2>=0;
pconstr(3) = S3>=0;
pconstr(4) = S4>=0;
pconstr(5) = S5>=0;
pconstr(6) = gamsym>=0;
pconstr(7) = 1-gamsym-1e-4>=0;
pconstr(8) = B>=0;
pconstr(9) = gamsym-1*B-S1*c(1)>=0;
pconstr(10) = B-1-S2*c(2)>=0;
pconstr(11) = -1*B+jacobian(B,x)*f+Mid-S3*c(4)>=0;

% Set object
obj = gamsym;

% Set options for sosopt
opts = sosoptions;
%opts.form = 'image';
opts.form = 'kernel';
opts.solver = 'mosek';
% opts.solveropts.epsilonStar = 1e-9;

% Call sosopt to find a barrier function
% [info,dopt,sossol] = sosopt(pconstr,x,opts);
[info,dopt,sossol] = sosopt(pconstr,x,obj,opts);

% Get barrier function
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Optimization is feasible.\n');
    Bs = subs(B,dopt);
    S1 = subs(S1,dopt);
    S2 = subs(S2,dopt);
    S3 = subs(S3,dopt);
    S4 = subs(S4,dopt);
    S5 = subs(S5,dopt);
    S = [S1,S2,S3,S4,S5];
else
    fprintf('\nSOS Optimization is not feasible. \n');    
    return
end

figure(2);clf;
xlim([-dom dom]); ylim([-dom dom]); hold on;   
[~,~]=pcontour(c1,0,domain,'m'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(c2,0,domain,'m'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(Bs,0,domain,'g'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(Bs,gam_0,domain,'r'); hold on;             % Plot the original Lyapunov sublevel set
Bdot = jacobian(Bs,x)*f;
[~,~]=pcontour(Bdot,0,domain,'b'); hold on;
axis(domain);

% end

% % Define SOS constraints
% % Constraint 1 :
% pconstr = -B+gam>=S2*L2;
% % Constraint 2: 
% pconstr(2) = B-1>=S3*L3;
% % Constraint 3: 
% pconstr(3) = jacobian(B,x)*f - B >=S4*L4;
% % Constraint 4: 
% pconstr(4) = B>=S5*L5;

% % Define SOS constraints
% % Constraint 1 :
% pconstr = -1*B+gam>=S2*c2;
% % Constraint 2: 
% pconstr(2) = B-1>=S3*c3;
% % Constraint 3: 
% pconstr(3) = B>=S4*c4;
% % Constraint 4: 
% pconstr(4) = jacobian(B,x)*f<=-S5*c5;
% % Constraint 5: 
% pconstr(5) = jacobian(B,x)*f>=S6*c3;

% % Define SOS constraints
% % Constraint 1 :
% pconstr = -1*B+gam>=S2*c2;
% % Constraint 2: 
% pconstr(2) = B-1>=S3*c3;
% % Constraint 3: 
% % pconstr(3) = B>=0;
% % Constraint 4: 
% Mid = 0.5'*0.5*sum(diag(jacobian(jacobian(B,x),x)))*0.5;
% pconstr(4) = jacobian(B,x)*f+Mid<=S5*c2;
% % Constraint 5: 
% % pconstr(5) = sum(diag(jacobian(jacobian(B,x),x)))>=0;

% % Define SOS constraints
% % Constraint 1 :
% pconstr = -1*B+gam>=S2*c2;
% % Constraint 2: 
% pconstr(2) = B-1>=S3*c3;
% % Constraint 3: 
% % pconstr(3) = B>=S4*c4;
% pconstr(3) = B>=0;
% % Constraint 4: 
% Mid = 0.5'*0.5*sum(diag(jacobian(jacobian(B,x),x)))*0.5;
% pconstr(4) = jacobian(B,x)*f+Mid<=S5*c2;
% % Constraint 5: 
% pconstr(5) = jacobian(B,x)*f+Mid<=S6*c3;

% % Verify the Gram Marix Decompositions in sossol
% z1 = sossol(1).z;
% Q1 = sossol(1).Q;
% e1 = (Bs-L1) - (z1'*Q1*z1);
% fprintf('\nMax magnitude coefficient of s(1)-z1''*Q1*z1 is:')
% disp(full(max(abs(e1.coefficient))))
% fprintf('Minimum eigenvalue of Q1 is:')
% disp(min(eig(Q1)));
% 
% z2 = sossol(2).z;
% Q2 = sossol(2).Q;
% e2 = (-jacobian(Bs,x)*f) - (z2'*Q2*z2);
% fprintf('\nMax magnitude coefficient of s(2)-z2''*Q2*z2 is:')
% disp(full(max(abs(e2.coefficient))))
% fprintf('Minimum eigenvalue of Q2 is:')
% disp(min(eig(Q2)));

