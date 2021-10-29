function [cc,kk,solu]=sos_function_v2(f,gg,k,k_l,V,C,dom,solL)

kk = 1;
domain = [-dom dom -dom dom];
pvar x1 x2 cc;
x = [x1;x2];
%     [u,uc] = polydecvar('u_w',monomials(x,0:k)); % L1 sos decision variables
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k_l)); % L1 sos decision variables
[L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:k_l)); % L1 sos decision variables
[L5,L5_Q] = sosdecvar('L5_w',monomials(x,0:k_l)); % L1 sos decision variables
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k)); % L1 sos decision variables

%     Vdot = jacobian(V, x1)*f(1) + jacobian(V, x2)*(f(2)+(x1+x2)*u);
Vdot = jacobian(V, x1)*f(1)+ jacobian(V, x2)*(f(2)+gg*u1);

%% Constraint:
pconstr_21 = L3 >= 0;
pconstr_22 = L4 >= 0;
pconstr_23 = L5 >= 0;
pconstr_1 = -Vdot-solL*(cc-V) >= 0;
pconstr_31 = -(cc-V)+C(1)*L3 >= 0;
pconstr_32 = -(cc-V)+C(2)*L4 >= 0;
pconstr_33 = -(cc-V)+C(3)*L5 >= 0;
pconstr_4 = cc >= 0;
% pconstr = [pconstr_21;pconstr_22;pconstr_23;pconstr_1;pconstr_31;pconstr_32;pconstr_33;pconstr_4];
pconstr = [pconstr_1;pconstr_4];

%% Set objection
obj = -cc;

%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);

figure(11);hold on;
% Create output
if info.feas
    cc = subs(cc,dopt);
    solu = subs(u1,dopt);
    [~,~]=pcontour(V,double(cc),domain,'g'); hold on;             % Plot the original Lyapunov sublevel set
    refreshdata; drawnow;
else
    kk = 0;
    cc = 0;
    solu = 0;
    fprintf('Barrier Certificate can not find.======\n');
    return;
end
end