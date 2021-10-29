function [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1_3D(f,k,solh,V,mm,gamma,gg,L_au,boundary_u)
pvar x1 x2 x3 htol;
x = [x1;x2;x3];
% Create corresponding decision variable
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:L_au/2)); % L1 sos decision variables
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:L_au/2)); % L1 sos decision variables
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k));
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k));
hdot = jacobian(solh, x1)*f(1) + jacobian(solh, x2)*(f(2)+gg(1)*u1)+jacobian(solh, x3)*(f(3)+gg(2)*u2);
Vdot = jacobian(V, x1)*f(1)+ jacobian(V, x2)*(f(2)+gg(1)*u1)+jacobian(V, x3)*(f(3)+gg(2)*u2);

%% Constrain :
sosconstr_1 = L1 >= 0;
sosconstr_2 = L2 >= 0;
sosconstr_3 = -Vdot-L1*solh >= 0;
sosconstr_4 = hdot+gamma*solh-L2*solh-htol >= 0;
sosconstr_5 = htol >=0;
sosconstr = [sosconstr_1;sosconstr_2;sosconstr_3;sosconstr_4;sosconstr_5];

input_con = [uc1;uc2];
for i=1:length(input_con)
    con = input_con(i);
    sosconstr = [sosconstr; con <= boundary_u];
    sosconstr = [sosconstr; con >= -boundary_u];
end

%% Set objection
obj = -htol;
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(sosconstr,x,obj,opts);
%% Create output
if info.feas
    SOL1 = subs(L1,dopt);
    SOL2 = subs(L2,dopt);
    SOLu1 = subs(u1,dopt);
    SOLu2 = subs(u2,dopt);
    kk = 1;
else
    kk = 0;
    SOL1 = 0;
    SOL2 = 0;
    SOLu1 = 0;
    SOLu2 = 0;
    fprintf('L1 and L2 can not find.====== ');
end
end