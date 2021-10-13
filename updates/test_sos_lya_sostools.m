% SOSDEMO2 --- Lyapunov Function Search 
% Section 3.2 of SOSTOOLS User's Manual
% 

clear; echo on;
% syms x1 x2 x3;
% vars = [x1; x2; x3];
syms x1 x2;       
vars = [x1;x2];

% Constructing the vector field dx/dt = f
% f = [-x1^3-x1*x3^2;
%     -x2-x1^2*x2;
%     -x3+3*x1^2*x3-3*x3/(x3^2+1)];
f = [x2-x1;
0.025359729472797157960711306552511*x1^6-0.01217009723261110433735733817297*x1^5-0.28561881441095026525735902480476*x1^4+0.0053102265416614026545590095169987*x1^3+x1^2*x2+0.17790534876627532116325861958709*x1^2-0.47421424943149480720390916606751*x1+0.083945383319193114801670674296474*x2^4+0.1356027624685912646995689101459*x2^3+0.071481552423721034239534333210031*x2^2-0.054409234652394423970012127256268*x2-0.000099625659358892892925041451235302
];

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vars);

% =============================================
% The Lyapunov function V(x): 
% [prog,V] = sospolyvar(prog,[x1^2; x2^2; x3^2],'wscoeff');
[prog,V] = sospolyvar(prog,[x1^2; x2^2],'wscoeff');

% =============================================
% Next, define SOSP constraints

% Constraint 1 : V(x) - (x1^2 + x2^2 + x3^2) >= 0
% prog = sosineq(prog,V-(x1^2+x2^2+x3^2));
prog = sosineq(prog,V-(x1^2+x2^2));

% Constraint 2: -dV/dx*(x3^2+1)*f >= 0
expr = -(diff(V,x1)*f(1)+diff(V,x2)*f(2))*(x2^2+1);
prog = sosineq(prog,expr);

% =============================================
% And call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

% =============================================
% Finally, get solution
SOLV = sosgetsol(prog,V)
echo off;