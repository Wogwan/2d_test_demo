function [solu1,solu2,solL,kk]=sos_function_v_3D(f,gg,k_u,k_L,V,cc,boundary_u)

kk = 1;
pvar x1 x2 x3;
x = [x1;x2;x3];
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k_u)); 
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k_u)); 
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k_L/2));
% Vdot = jacobian(V, x1)*(x2-x1)+ jacobian(V, x2)*(0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608+gg*u);
Vdot = jacobian(V, x1)*f(1)+ jacobian(V, x2)*(f(2)+gg(1)*u1)+jacobian(V, x3)*(f(3)+gg(2)*u2);

%% Constraint:
pconstr_1 = L3 >= 0;
pconstr_2 = -Vdot-L3*(cc-V) >= 0;
pconstr = [pconstr_1; pconstr_2];
% pconstr_3 = sum([uc1;uc2])/length([uc1;uc2]) <= 100;
input_con = [uc1;uc2];
for i=1:length(input_con) 
        con = input_con(i);
        pconstr = [pconstr; con <= boundary_u];
        pconstr = [pconstr; con >= -boundary_u];
end

%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
%     [info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
[info,dopt] = sosopt(pconstr,x,opts);

% Create output
if info.feas
    solu1 = subs(u1,dopt);
    solu2 = subs(u2,dopt);
    solL = subs(L3,dopt);
else
    kk = 0;
    solu1 = 0;
    solu2 = 0;
    solL = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end
end