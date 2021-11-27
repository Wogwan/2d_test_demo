function [solu1,solu2,solu3,solL,kk]=sos_function_v_3D(f,gg,k_u,k_L,V,cc)
kk = 1;
pvar x1 x2 x3;
x = [x1;x2;x3];
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k_u));
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k_u));
[u3,uc3] = polydecvar('u_w3',monomials(x,0:k_u));
%%
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k_L/2));
%%
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2)+jacobian(V, x3)*(f(3)+gg(3)*u3);
%% Constraint:
pconstr_1 = L3 >= 0;
pconstr_2 = -Vdot-L3*(cc-V) >= 0;
pconstr = [pconstr_1; pconstr_2];
%% Input Limits
% input_con = [uc1;uc2];
% for i=1:length(input_con)
%         con = input_con(i);
%         pconstr = [pconstr; con <= boundary_u];
%         pconstr = [pconstr; con >= -boundary_u];
% end
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,opts);
%% Create output
if info.feas
    solu1 = subs(u1,dopt);
    solu2 = subs(u2,dopt);
    solu3 = subs(u3,dopt);
    solL = subs(L3,dopt);
else
    kk = 0;
    solu1 = 0;
    solu2 = 0;
    solu3 = 0;
    solL = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end
end