function [Bs, Gam] = barrier_gam(k,f,B,S,c)

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
pvar x1 x2 gam;
x = [x1;x2];
% [B,Bc] = polydecvar('b',monomials(x,0:k));

% [S2,S2_Q] = sosdecvar('S2_w',monomials(x,0:k)); % S2 sos decision variables
% [S3,S3_Q] = sosdecvar('S3_w',monomials(x,0:k)); % S3 sos decision variables
% [S4,S4_Q] = sosdecvar('S4_w',monomials(x,0:k)); % S4 sos decision variables
% [S5,S5_Q] = sosdecvar('S5_w',monomials(x,0:k)); % S5 sos decision variables
% [S6,S6_Q] = sosdecvar('S6_w',monomials(x,0:k)); % S6 sos decision variables

% Define SOS constraints
% Constraint 1 :
pconstr = -1*B+gam>=S(1)*c(1);
% Constraint 2: 
pconstr(2) = B-1>=S(2)*c(2);
% Constraint 3: 
pconstr(3) = B>=S(3)*c(3);
% pconstr(3) = B>=0;
% Constraint 4: 
pconstr(4) = 1-B>=S(4)*c(4);
% Constraint 5: 
Mid = 0.5'*0.5*sum(diag(jacobian(jacobian(B,x),x)))*0.5;
pconstr(5) = jacobian(B,x)*f+Mid<=S(5)*c(1);
% Constraint 6: 
% pconstr(6) = jacobian(B,x)*f+Mid<=S6*c3;

% Set object
obj = -gam;

% Set options for sosopt
opts = sosoptions;
%opts.form = 'image';
opts.form = 'kernel';
opts.solver = 'mosek';
% opts.solveropts.epsilonStar = 1e-9;

% Call sosopt to find a barrier function
[info,dopt,sossol] = sosopt(pconstr,x,obj,opts);
% [info,dopt,sossol] = sosopt(pconstr,x,opts);

% Get barrier function
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Gam Optimization is feasible.\n');
    Gam = subs(gam,dopt);
%     Bs = subs(B,dopt);
else
    fprintf('\nSOS Gam Optimization is not feasible. \n');    
    return
end

end