function [cc,kk,solu]=sos_function_v2_3D(f,gg,k_u,k_L,V,C,dom,solL,boundary_u)

kk = 1;
domain = [-dom dom -dom dom -dom dom];
pvar x1 x2 x3 cc;
x = [x1;x2;x3];
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k_L/2)); 
[L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:k_L/2)); 
[L5,L5_Q] = sosdecvar('L5_w',monomials(x,0:k_L/2)); 
[L6,L6_Q] = sosdecvar('L6_w',monomials(x,0:k_L/2)); 
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k_u/2)); 
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k_u/2)); 
%     Vdot = jacobian(V, x1)*f(1) + jacobian(V, x2)*(f(2)+(x1+x2)*u);
Vdot = jacobian(V, x1)*f(1)+ jacobian(V, x2)*(f(2)+gg(1)*u1)+ jacobian(V, x3)*(f(3)+gg(2)*u2);

%% Constraint:
pconstr_21 = L3 >= 0;
pconstr_22 = L4 >= 0;
pconstr_23 = L5 >= 0;
pconstr_24 = L6 >= 0;
pconstr_1 = -Vdot-solL*(cc-V) >= 0;
pconstr_31 = -(cc-V)+C(1)*L3 >= 0;
pconstr_32 = -(cc-V)+C(2)*L4 >= 0;
pconstr_33 = -(cc-V)+C(3)*L5 >= 0;
pconstr_34 = -(cc-V)+C(4)*L6 >= 0;
pconstr_4 = cc >= 0;
pconstr = [pconstr_21;pconstr_22;pconstr_23;pconstr_24;pconstr_1;pconstr_31;pconstr_32;pconstr_33;pconstr_34;pconstr_4];
% pconstr = [pconstr_21;pconstr_22;pconstr_24;pconstr_1;pconstr_31;pconstr_32;pconstr_34;pconstr_4];
% pconstr = [pconstr_1;pconstr_4];

input_con = [uc1;uc2];
for i=1:length(input_con) 
        con = input_con(i);
        pconstr = [pconstr; con <= boundary_u];
        pconstr = [pconstr; con >= -boundary_u];
end

%% Set objection
obj = -cc;

%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,obj,opts);

figure(11);hold on;
% Create output
if info.feas
    cc = subs(cc,dopt);
    solu = subs(u1,dopt);
    inV = patch(pcontour3(V,double(cc),domain,'k')); 
    set(inV, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',0.7 ); hold on;
    refreshdata; drawnow;
else
    kk = 0;
    cc = 0;
    solu = 0;
    fprintf('Barrier Certificate can not find.======\n');
    return;
end
end