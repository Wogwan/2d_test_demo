function [cc,kk,solu]=sos_function_v2_3D(f,gg,k_u,k_L,V,C,dom,solL,figure_id)

kk = 1;
domain = [-dom dom -dom dom -dom dom];
pvar x1 x2 x3 cc;
x = [x1;x2;x3];
%%
% [L3,L3_Q] = polydecvar('L1_w',monomials(x,0:k_L));
% [L4,L4_Q] = polydecvar('L1_w',monomials(x,0:k_L)); 
% [L5,L5_Q] = polydecvar('L1_w',monomials(x,0:k_L)); 
% [L6,L6_Q] = polydecvar('L1_w',monomials(x,0:k_L));
%%
[L3,L3_Q] = sosdecvar('L1_w',monomials(x,0:k_L/2));
[L4,L4_Q] = sosdecvar('L1_w',monomials(x,0:k_L/2)); 
[L5,L5_Q] = sosdecvar('L1_w',monomials(x,0:k_L/2)); 
[L6,L6_Q] = sosdecvar('L1_w',monomials(x,0:k_L/2));
%%
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k_u)); 
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k_u)); 
[u3,uc3] = polydecvar('u_w3',monomials(x,0:k_u)); 
%     Vdot = jacobian(V, x1)*f(1) + jacobian(V, x2)*(f(2)+(x1+x2)*u);
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2)+jacobian(V, x3)*(f(3)+gg(3)*u3);

%% Constraint:
pconstr_21 = L3 >= 0;
pconstr_22 = L4 >= 0;
pconstr_23 = L5 >= 0;
pconstr_24 = L6 >= 0;
pconstr_1 = -Vdot-solL*(cc-V) >= 0;
pconstr_31 = -(cc-V)+C(1)*L3 >= 0;
pconstr_32 = -(cc-V)+C(2)*L4 >= 0;
pconstr_33 = -(cc-V)+C(3)*L5 >= 0;
% pconstr_34 = -(cc-V)+C(4)*L6 >= 0;
pconstr_4 = cc >= 0;
% pconstr = [pconstr_21;pconstr_22;pconstr_23;pconstr_24;pconstr_1;pconstr_31;pconstr_32;pconstr_33;pconstr_34;pconstr_4];
pconstr = [pconstr_21;pconstr_22;pconstr_23;pconstr_1;pconstr_31;pconstr_32;pconstr_33;pconstr_4];

%% Input limits
% input_con = [uc1;uc2];
% for i=1:length(input_con) 
%         con = input_con(i);
%         pconstr = [pconstr; con <= boundary_u];
%         pconstr = [pconstr; con >= -boundary_u];
% end
%% Set objection
obj = -cc;
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,obj,opts);
figure(figure_id);hold on;
% Create output
if info.feas
    cc = subs(cc,dopt);
    solu1 = subs(u1,dopt);
    solu2 = subs(u2,dopt);
    solu3 = subs(u3,dopt);
    solu = [solu1;solu2;solu3];
    inV = patch(pcontour3(V,double(cc),domain,'k')); 
    set(inV, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',0.7 ); hold on;
    refreshdata; drawnow;
else
    kk = 0;
    cc = 0;
    solu = [0;0;0];
    fprintf('Barrier Certificate can not find.======\n');
    return;
end
end