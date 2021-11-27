function [solL,kk]=sos_optimal_v2_3D(f,gg,k,k_l,V,cc)

kk = 1;
pvar x1 x2 x3;
x = [x1;x2;x3];
%%
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k));
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k));
[u3,uc3] = polydecvar('u_w3',monomials(x,0:k));
%%
[L ,L_Q] = polydecvar('L_w',monomials(x,0:k_l));
%%
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+jacobian(V, x2)*(f(2)+gg(2)*u2)+jacobian(V, x3)*(f(3)+gg(3)*u3);
%% Constraint:
pconstr_1 = L >= 0;
pconstr_2 = -Vdot-L*(cc-V) >= 0;
pconstr = [pconstr_1; pconstr_2];
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,opts);
%% Create output
if info.feas
    solL = subs(L,dopt);
else
    kk = 0;
    solL = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end
end