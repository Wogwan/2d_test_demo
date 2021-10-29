function [V,kk]=sos_optimal_v(f,gg,k,B,u,C)

kk = 1;
pvar x1 x2 Vtol;
%%
% f = [x2-x1
%     0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608
%     ];
% gg = 1;
% k = 2;
% dom = 8;
% B = -20920.02735971049*x1^2-12944.00205666696*x1*x2-21004.89487163851*x2^2+725.8021517955463*x1-34941.78371754049*x2+132960.4633676535;
% u = -1.423136761279898*x1^2+0.1389802653220975*x1*x2+0.1251182246445742*x2^2-1.103317535267957*x1-1.214738774413024*x2-0.02619321979348357;

%%
x = [x1;x2];
[V,vc] = polydecvar('v_w',monomials(x,0:k)); % L1 sos decision variables
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:k)); % L1 sos decision variables
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:k)); % L1 sos decision variables
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k)); % L1 sos decision variables
[L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:k)); % L1 sos decision variables
[L5,L5_Q] = sosdecvar('L5_w',monomials(x,0:k)); % L1 sos decision variables

% Vdot = jacobian(V, x1)*(x2-x1)+ jacobian(V, x2)*(0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608+gg*u);
Vdot = jacobian(V, x1)*f(1)+ jacobian(V, x2)*(f(2)+gg*u);

%% Constraint:
% pcr_1 = L1 >= 0;
% pcr_2 = L2 >= 0;
% pconstr_1 = -V + L1*B >= 0;
% pconstr_2 = -Vdot-L2*B-Vtol >= 0;
% pconstr_3 = Vtol >= 0;
% pconstr = [pcr_1; pcr_2; pconstr_1; pconstr_2; pconstr_3];

% pcr_1 = L1 >= 0;
% pcr_2 = L2 >= 0;
% pconstr_1 = -V + L1*B >= 0;
% pconstr_2 = Vdot-L2*B+Vtol >= 0;
% pconstr_3 = Vtol >= 0;
% pconstr = [pcr_1; pcr_2; pconstr_1; pconstr_2; pconstr_3];

pcr_1 = L1 >= 0;
pcr_2 = L2 >= 0;
pcr_6 = L3 >= 0;
pcr_7 = L4 >= 0;
pcr_8 = L5 >= 0;
pcr_31 = -V+C(1)*L3 >= 0;
pcr_32 = -V+C(2)*L4 >= 0;
pcr_33 = -V+C(3)*L5 >= 0;
% pcr_31 = V+C(1)*L3 >= 0;
% pcr_32 = V+C(2)*L4 >= 0;
% pcr_33 = V+C(3)*L5 >= 0;
pconstr_1 = V-L1*B >= 0;
pconstr_2 = Vdot+L2*B-Vtol >= 0;
pconstr_3 = Vtol >= 0;
% pconstr = [pcr_1; pcr_2; pcr_6; pcr_7; pcr_8; pcr_31; pcr_32; pcr_33; pconstr_1; pconstr_2; pconstr_3];
pconstr = [pconstr_1; pconstr_2; pconstr_3];

%% Set objection
obj = Vtol;

%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
%     [info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
[info,dopt] = sosopt(pconstr,x,obj,opts);

% Create output
if info.feas
    V = subs(V,dopt);
else
    kk = 0;
    V  = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end
end