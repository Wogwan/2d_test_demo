%---------------------------------------------------------------------
% sosoptdemo2
%
% Demonstration of SOSOPT function for solving sos optimizations.
% This example uses SOSOPT to find a Lyapunov function for a three
% state rational system. 
% See SOSOPT help for more details on the function syntax.
%
% Reference: 
% S. Prajna, A. Papachristodoulou, P. Seiler, and P.A. Parrilo,
% SOSTOOLS: SOSDEMO2, Section 3.2 of SOSTOOLS Version 2.00
% User's Guide, 2004.
%---------------------------------------------------------------------

format long

% Define polynomial variables
% pvar x1 x2 x3;       
% x = [x1;x2;x3];
pvar x1 x2;       
x = [x1;x2];


% Create vector field: dx/dt = f(x)/(x3^2+1)
% f = [(-x1^3-x1*x3^2)*(x3^2+1);
%     (-x2-x1^2*x2)*(x3^2+1);
%     (-x3+3*x1^2*x3)*(x3^2+1)-3*x3];

f = [x2-x1;
0.025359729472797157960711306552511*x1^6-0.01217009723261110433735733817297*x1^5-0.28561881441095026525735902480476*x1^4+0.0053102265416614026545590095169987*x1^3+x1^2*x2+0.17790534876627532116325861958709*x1^2-0.47421424943149480720390916606751*x1+0.083945383319193114801670674296474*x2^4+0.1356027624685912646995689101459*x2^3+0.071481552423721034239534333210031*x2^2-0.054409234652394423970012127256268*x2-0.000099625659358892892925041451235302
];

% Create Lyapunov function V(x) with unspecified coefficients
% zV = [x1^2; x2^2; x3^2];
% [V,cV] = polydecvar('c',zV);
% zV = [x1; x2; x1*x2; x1^2; x2^2];
zV = [x1^2*x2^2; x1^4; x2^4; x1^2; x2^2];
[V,cV] = polydecvar('c',zV);

% Define SOS constraints
% Constraint 1 : V(x) >= (x1^2 + x2^2 + x3^2) 
L1 = x'*x;
pconstr = V>=0;

% Constraint 2: dV/dt = dV/dx*(x3^2+1)*f <= 0
pconstr(2) = jacobian(V,x)*f<=0;

% Set options for sosopt
opts = sosoptions;
% opts.form = 'image';
opts.form = 'kernel';
opts.solver = 'mosek';
%opts.solver = 'sdpt3';
%opts.solver = 'csdp';
%opts.solver = 'sdplr';
%opts.solver = 'sdpam'; opts.solveropts.epsilonStar = 1e-9;

% Call sosopt to find a Lyapunov function
[info,dopt,sossol] = sosopt(pconstr,x,opts);

% Get Lyapunov function
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Optimization is feasible.\n');
    Vs = subs(V,dopt)
else
    fprintf('\nSOS Optimization is not feasible. \n');    
    return
end

% Verify the Gram Marix Decompositions in sossol
z1 = sossol(1).z;
Q1 = sossol(1).Q;
e1 = (Vs-L1) - (z1'*Q1*z1);
fprintf('\nMax magnitude coefficient of s(1)-z1''*Q1*z1 is:')
disp(full(max(abs(e1.coefficient))))
fprintf('Minimum eigenvalue of Q1 is:')
disp(min(eig(Q1)));

z2 = sossol(2).z;
Q2 = sossol(2).Q;
e2 = (-jacobian(Vs,x)*f) - (z2'*Q2*z2);
fprintf('\nMax magnitude coefficient of s(2)-z2''*Q2*z2 is:')
disp(full(max(abs(e2.coefficient))))
fprintf('Minimum eigenvalue of Q2 is:')
disp(min(eig(Q2)));


