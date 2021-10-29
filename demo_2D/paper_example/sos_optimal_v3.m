function [V,kk]=sos_optimal_v3(f,gg,k,L_au,B,u)

kk = 1;
pvar x1 x2 Vtol;
%%
f = [x2-x1
    0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608
    ];
gg = 1;
k = 2;
dom = 8;
solh = -20920.02735971049*x1^2-12944.00205666696*x1*x2-21004.89487163851*x2^2+725.8021517955463*x1-34941.78371754049*x2+132960.4633676535;
u = -1.423136761279898*x1^2+0.1389802653220975*x1*x2+0.1251182246445742*x2^2-1.103317535267957*x1-1.214738774413024*x2-0.02619321979348357;
L_au = 2;

%%
x = [x1;x2];
% Create corresponding decision variable
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:L_au)); % L1 sos decision variables
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:L_au)); % L1 sos decision variables
[u,u1_Q] = polydecvar('u1_w',monomials(x,0:k)); % u1 sos decision variables
[V,vc] = polydecvar('v_w',monomials(x,0:k)); % L1 sos decision variables

hdot = jacobian(solh, x1)*f(1) + jacobian(solh, x2)*(f(2)+gg*u);
Vdot = jacobian(V, x1)*f(1) + jacobian(V, x2)*(f(2)+gg*u);
%%
% Constrain :
sosconstr_1 = L1 >= 0;
sosconstr_2 = L2 >= 0;
sosconstr_3 = -Vdot-L1*solh >= 0;
sosconstr_4 = V - L2*solh >= Vtol;
sosconstr_5 = Vtol >=0;
sosconstr = [sosconstr_1;sosconstr_2;sosconstr_3;sosconstr_4;sosconstr_5];
%% Set objection
obj = Vtol;
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(sosconstr,x,obj,opts);
%% Create output
if info.feas
    SOL1 = subs(L1,dopt);
    SOL2 = subs(L2,dopt);
    SOLu = subs(u,dopt);
    kk = 1;
else
    SOL1 = 0;
    SOL2 = 0;
    SOLu = 0;
    kk = 0;
    fprintf('L1 and L2 can not find.====== ');
end
end