function [V, kk] = sos_optimal_V1(f,gg,B,u1,u2,l_au,l_us,V_degree,C,gamma)

pvar x1 x2 Vtol2;
x = [x1;x2];
%%
[V,vc] = polydecvar('v_w',monomials(x,0:V_degree));
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:l_au));
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:l_au)); 
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:l_us)); 
[L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:l_us)); 
%%
% hdot = jacobian(B, x1)*(f(1)+gg(1)*u1)+ jacobian(B, x2)*(f(2)+gg(2)*u2);
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
%% Constraint
pcr_11 = L1 >= 0;
pcr_12 = L2 >= 0;
pcr_13 = L3 >= 0;
pcr_14 = L4 >= 0;
pcr_21 = V-C(1)*L3 >= 0;
pcr_22 = V-C(2)*L4 >= 0;
pconstr_1 = V-L1*B >= 0;
pconstr_2 = -Vdot-L2*B+gamma*B-Vtol2 >= 0;
pconstr_3 = Vtol2 >= 0;
pconstr = [pcr_11;pcr_12;pcr_13;pcr_14;pconstr_1;pconstr_2;pconstr_3;pcr_21;pcr_22];
%% Set objection
obj = Vtol2;

% pconstr_2 = -Vdot-L2*B+gamma*B >= 0;
% pconstr_3 = Vtol2 >= 0;
% pconstr = [pcr_11;pcr_12;pcr_13;pcr_14;pconstr_1;pconstr_2;pcr_21;pcr_22];

%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,obj,opts);
% [info,dopt] = sosopt(pconstr,x,opts);

% Create output
if info.feas
    kk = 1;
    V = subs(V,dopt)
else
    kk = 0;
    V  = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end
end