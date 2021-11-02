function [cc,kk,solu1,solu2]=sos_optimal_v3(f,gg,k,k_l,V,C,dom,solL,ccc,figure_id)

kk = 1;
domain = [-dom dom -dom dom];
pvar x1 x2 cc;
x = [x1;x2];
%     [u,uc] = polydecvar('u_w',monomials(x,0:k)); % L1 sos decision variables
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:k_l/2)); % L1 sos decision variables
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:k_l/2)); % L1 sos decision variables
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k_l/2)); % L1 sos decision variables
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k)); % L1 sos decision variables
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k)); % L1 sos decision variables

% Vdot = jacobian(V, x1)*f(1) + jacobian(V, x2)*(f(2)+gg(2)*u1);
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);

%% Constraint:
pconstr_21 = L1 >= 0;
pconstr_22 = L2 >= 0;
% pconstr_23 = L3 >= 0;
pconstr_1 = -Vdot-solL*(cc-V) >= 0;
pconstr_31 = -(cc-V)+C(1)*L1 >= 0;
pconstr_32 = -(cc-V)+C(2)*L2 >= 0;
% pconstr_33 = -(cc-V)+C(3)*L3 >= 0;
pconstr_4 = cc >= ccc;
pconstr = [pconstr_21;pconstr_22;pconstr_1;pconstr_31;pconstr_32;pconstr_4];
% pconstr = [pconstr_1;pconstr_4];

%% Set objection
obj = -cc;

%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
%%
figure(figure_id);hold on;
if info.feas
    cc = subs(cc,dopt);
    solu1 = subs(u1,dopt);
    solu2 = subs(u2,dopt);
    [~,~]=pcontour(V,double(cc),domain,'r'); hold on;             % Plot the original Lyapunov sublevel set
    refreshdata; drawnow;
else
    kk = 0;
    cc = 0;
    solu1 = 0;
    solu2 = 0;
    fprintf('Suitable sublevel set can not find.======\n');
    return;
end
end