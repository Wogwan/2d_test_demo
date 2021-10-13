function [Bs, S] = barrier_test(k,f,gam_0,c)

%%
% k = 4;
% dom = 4;
% domain = [-dom dom -dom dom];
% gam_0 = 0.1;
% c2 = gam_0 -(x1^2+x2^2); % Safe region qualified polynomial X_0
% c3 = x1^2+x2^2-4;    % Unsafe region : X_u
% c4 = x1^2+x2^2;       % For all [x x] : X
% c5 = c4-c3;              % X \ X_u
% f = [x2; 
%     (1-x1^2)*x2-x1];

%%
% Define polynomial variables
pvar x1 x2;       
x = [x1;x2];

% [B,Bc] = polydecvar('b',monomials(x,0:k));
[B,Bc] = sosdecvar('b',monomials(x,0:k));

[S2,S2_Q] = sosdecvar('S2_w',monomials(x,0:k)); % S2 sos decision variables
[S3,S3_Q] = sosdecvar('S3_w',monomials(x,0:k)); % S3 sos decision variables
[S4,S4_Q] = sosdecvar('S4_w',monomials(x,0:k)); % S4 sos decision variables
[S5,S5_Q] = sosdecvar('S5_w',monomials(x,0:k)); % S5 sos decision variables
[S6,S6_Q] = sosdecvar('S6_w',monomials(x,0:k)); % S6 sos decision variables

% Define SOS constraints
% Constraint 1 :
pconstr = -1*B+gam_0>=S2*c(1);
% Constraint 2: 
pconstr(2) = B-1>=S3*c(2);
% Constraint 3: 
pconstr(3) = B>=S4*c(3);
% pconstr(3) = B>=0;
% Constraint 4: 
pconstr(4) = 1-B>=S5*c(4);
% Constraint 5: 
Mid = 0.5'*0.5*sum(diag(jacobian(jacobian(B,x),x)))*0.5;
pconstr(5) = jacobian(B,x)*f+Mid<=S6*c(1);
% Constraint 6: 
% pconstr(6) = jacobian(B,x)*f+Mid<=S6*c3;

% Set options for sosopt
opts = sosoptions;
%opts.form = 'image';
opts.form = 'kernel';
opts.solver = 'mosek';
% opts.solveropts.epsilonStar = 1e-9;

% Call sosopt to find a barrier function
[info,dopt,sossol] = sosopt(pconstr,x,opts);

% Get barrier function
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Optimization is feasible.\n');
    Bs = subs(B,dopt);
    S2 = subs(S2,dopt);
    S3 = subs(S3,dopt);
    S4 = subs(S4,dopt);
    S5 = subs(S5,dopt);
    S6 = subs(S6,dopt);
    S = [S2,S3,S4,S5,S6];
else
    fprintf('\nSOS Optimization is not feasible. \n');    
    return
end

% figure(2);clf;
% xlim([-dom dom]); ylim([-dom dom]); hold on;   
% [~,~]=pcontour(c2,0,domain,'m'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(c3,0,domain,'m'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(Bs,0,domain,'g'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(Bs,gam_0,domain,'r'); hold on;             % Plot the original Lyapunov sublevel set
% Bdot = jacobian(Bs,x)*f;
% [~,~]=pcontour(Bdot,0,domain,'b'); hold on;
% axis(domain);

end

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

