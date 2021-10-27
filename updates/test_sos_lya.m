% Define polynomial variables
format long
pvar x1 x2;
x = [x1;x2];
% Create vector field: dx/dt = f(x)/(x3^2+1)
f = [x2-x1
    0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608
    ];

% Create Lyapunov function V(x) with unspecified coefficients
zV = [x1^2*x2^2; x1^4; x2^4; x1^2; x2^2;x1^2*x2^2];
[V,cV] = polydecvar('c',zV);

% Define SOS constraints
% Constraint 1 : V(x) >= (x1^2 + x2^2 + x3^2)
L1 = x'*x;
pconstr = V>=0;

% Constraint 2: dV/dt = dV/dx*(x3^2+1)*f <= 0
pconstr(2) = jacobian(V,x)*f<=0;

% Set options for sosopt
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
% opts.solveropts.epsilonStar = 1e-9;

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


