function [solu,solL,kk]=sos_function_v(f,gg,k,V,cc,dom)

kk = 1;
pvar x1 x2 u;
x = [x1;x2];
[u,uc] = polydecvar('u_w',monomials(x,0:k)); % L1 sos decision variables
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k/2)); % L1 sos decision variables
% Vdot = jacobian(V, x1)*(x2-x1)+ jacobian(V, x2)*(0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608+gg*u);
    Vdot = jacobian(V, x1)*f(1)+ jacobian(V, x2)*(f(2)+gg*u);

%% Constraint:
pconstr_1 = L3 >= 0;
pconstr_2 = -Vdot-L3*(cc-V) >= 0;
pconstr = [pconstr_1; pconstr_2];

%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
%     [info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
[info,dopt] = sosopt(pconstr,x,opts);

% Create output
if info.feas
    solu = subs(u,dopt);
    solL = subs(L3,dopt);
else
    kk = 0;
    solu  = 0;
    solL = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end
end