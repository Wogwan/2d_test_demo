function k_iter = auto_test_demo_2d(sym,lp1,lp2,lp3,lp4,lp5,plot)
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
k_iter = 0;
%% Initial System information
% C1 = (x1+8)^2+(x2-0)^2-4;
% C2 = (x1-8)^2+(x2+0)^2-4;
% C3 = (x1-0)^2+(x2-8)^2-4;
% C4 = (x1-0)^2+(x2+8)^2-4;
% sym.us = [C1;C2;C3;C4];
% % sym.us = [C1;C2;C3];
% sym.f = [x2-x1;
%     0.46382156330341729589940137480476*x1^4-0.21219821363526560116570580989732*x1^3+0.35000677410250802565647474138031*x1^2+x1*x2-0.27258976376810471290805063896793*x1+0.15407267311901634565529661813343*x2^4+0.94268273772394006737584959410015*x2^3+4.1277331008261466394060335005634*x2^2-1.7098389065959235244562819389103*x2+0.011510115809394316777058975276304
%     ];
% sym.gg = [1; 1];
% sym.sys = [sym.f(1)+sym.gg(1)*u1; sym.f(2)+sym.gg(2)*u2];
% sym.V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
% % sym.V0 = 4*x1^4+2*x2^4+2*x1^2*x2^2+4*x1^2+2*x2^2+1*x1*x2;
% % sym.V0 = x1^2+x1*x2+x2^2;

%%
N_Lya = [];
lp1.solU = []; lp1.v_c = [];
TRACE = []; Barrier = []; Control = [];
TRACE_lp5 = []; Barrier_lp5 = []; Control_lp5 = [];

%% Plot basic figure
lp3.dom_2 = 100; domain_2 = [-lp3.dom_2 lp3.dom_2 -lp3.dom_2 lp3.dom_2];
dom = 15; plot.domain = [-dom dom -dom dom];
% plot.figure_id = 111;
figure(plot.figure_id);clf;hold on;
[~,~]=pcontour(sym.us(1),0,plot.domain,'r');
[~,~]=pcontour(sym.us(2),0,plot.domain,'r');
if length(sym.us) == 3
    [~,~]=pcontour(sym.us(3),0,plot.domain,'r');
elseif length(sym.us) == 4
    [~,~]=pcontour(sym.us(3),0,plot.domain,'r');
    [~,~]=pcontour(sym.us(4),0,plot.domain,'r');
else
    fprintf('40 The constraint number does not match.======\n');
end

%% Hyperparameter of lp1 @ CLF sublevel set
% lp1.v_k_u = 2;
% lp1.v_k_l = 4;
%%
lp1.iter = 1;
lp1.feas = 1;
lp1.epsi = 1e-6;
lp1.c0 = 1;
lp1.cc = 1.1;
[~,~]=pcontour(sym.V0,lp1.c0,plot.domain,'g');
%% Start to compute the sublevel set of a Lyapunov function
while double(lp1.cc)-double(lp1.c0) >= lp1.epsi
    lp1.iter = lp1.iter + 1;
    if lp1.iter ~= 1
        lp1.c0 = lp1.cc;
    end
    [solu1,solu2,solL,lp1.feas]= sos_function_v(sym.f,sym.gg,lp1.v_k_u,lp1.v_k_l,sym.V0,lp1.c0);
    if lp1.feas == 0
        break
    end
    [lp1.cc,lp1.feas,solu1,solu2] = sos_function_v2(sym.f,sym.gg,lp1.v_k_u,lp1.v_k_l,sym.V0,sym.us,dom,solL,plot.figure_id);
    double(lp1.cc)
    lp1.v_c = [lp1.v_c; double(lp1.cc)];
    lp1.solU = [lp1.solU;[solu1,solu2]];
    if lp1.feas == 0
        break
    end
end
hold on;
[~,~]=pcontour(sym.V0,max(double(lp1.v_c)),plot.domain,'b');

%% Plot
figure(plot.figure_id+1);clf;hold on;
[~,~]=pcontour(sym.us(1),0,plot.domain,'r');
[~,~]=pcontour(sym.us(2),0,plot.domain,'r');
if length(sym.us) == 3
    [~,~]=pcontour(sym.us(3),0,plot.domain,'r');
elseif length(sym.us) == 4
    [~,~]=pcontour(sym.us(3),0,plot.domain,'r');
    [~,~]=pcontour(sym.us(4),0,plot.domain,'r');
else
    fprintf('84 The constraint number does not match.======\n');
end
%% Hyperparameters of the lp2 @ CBF find controller
% lp2.u = 4;
% lp2.h = 4;
% lp2.us = 8;
% lp2.au = 8;
% lp2.gamma = 0;
%%
lp2.c_b = max(double(lp1.v_c));
lp2.solh = lp2.c_b - sym.V0;
lp2.epsi2 = 0.196;
lp2.feas = 1; j = 0;
%% Start to compute the control barrier function
record_Q = [1000];
trace_Q = 1000;
while abs(double(trace_Q)-double(record_Q(end)))>= lp2.epsi2
    j = j+1
    record_Q = [record_Q; trace_Q];
    record_Q
    [SOLu1,SOLu2,SOL1,SOL2,lp2.feas] = sos_function_1(sym.f,lp2.u,lp2.au,lp2.solh,sym.V0,lp2.gamma,sym.gg);
    if lp2.feas == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2]];
    [solh,trace_Q,lp2.feas]=sos_function_2(j,sym.f,lp2.h,SOLu1,SOLu2,SOL1,SOL2,lp2.gamma,sym.V0,sym.us,dom,sym.gg,lp2.us,plot.figure_id+1);
    if lp2.feas == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
end
figure(plot.figure_id);hold on;
if ~isempty(Barrier)
    [~,~]=pcontour(Barrier(end),0,plot.domain,'m');
else
    fprintf('120 Barrier function can not find.======\n');
end
%% Hyperparameters of the SOSP @ CBF -> V
% lp3.V_us = 6;
% lp3.V_au = 8;
% lp3.V_degree = 4;
% lp3.gamma = 0;
%%
for 1
    kk = 1; OO = 0;
    if ~isempty(Control)
        u1 = Control(end,1);
        u2 = Control(end,2);
    else
        u1 = 1;
        u2 = 1;
        B = 0;
        fprintf('137 Barrier function can not find.======\n');
    end
    if ~isempty(Barrier)
        B = Barrier(end);
    else
        B = lp2.c_b - sym.V0;
    end
    [lp3.V, lp3.feas] = sos_optimal_V1(sym.f,sym.gg,B,u1,u2,lp3.V_au,lp3.V_us,lp3.V_degree,sym.us,lp3.gamma);
    if lp3.feas == 0
        fprintf('146 Suitable Lyapunov function can not find.======\n');
    end
    N_Lya = [N_Lya;lp3.V];
    figure(plot.figure_id+2);clf;hold on;
    cc_p = [1,2,3];
    for i = 0:16
        k_u_V = ['r','g','b','m','c','k','y'];
        if mod(i,7) == 0
            [~,~]=pcontour(lp3.V,cc_p(1)*i,domain_2,k_u_V(7)); hold on;
            [~,~]=pcontour(lp3.V,cc_p(2)*i,domain_2,k_u_V(7)); hold on;
            [~,~]=pcontour(lp3.V,cc_p(3)*i,domain_2,k_u_V(7)); hold on;
        else
            [~,~]=pcontour(lp3.V,cc_p(1)*i,domain_2,k_u_V(mod(i,7))); hold on;
            [~,~]=pcontour(lp3.V,cc_p(2)*i,domain_2,k_u_V(mod(i,7))); hold on;
            [~,~]=pcontour(lp3.V,cc_p(3)*i,domain_2,k_u_V(mod(i,7))); hold on;
        end
    end
    lp3.Vdot = jacobian(lp3.V, x1)*(sym.f(1)+sym.gg(1)*u1)+ jacobian(lp3.V, x2)*(sym.f(2)+sym.gg(2)*u2);
    for i = 0:16
        k_u_V = ['r','g','b','m','c','k','y'];
        if mod(i,7) == 0
            [~,~]=pcontour(lp3.Vdot,0,domain_2,k_u_V(7)); hold on;             % Plot the original Lyapunov sublevel set
        else
            [~,~]=pcontour(lp3.Vdot,0,domain_2,k_u_V(mod(i,7))); hold on;             % Plot the original Lyapunov sublevel set
        end
    end
    %% Hyperparameters of the Optimal Sublevel set @ CLF LP4
    %     lp4.u = 4;
    %     lp4.au = 8;
    %%
    lp4.u1 = u1;
    lp4.u2 = u2;
    lp4.feas = 1;
    [a1,b1] = coeffs(p2s(lp3.V));
    ccc = double(vpa(a1(end)));
    lp4.V = lp3.V;
    C0 = ccc; CC = ccc+0.1;
    figure(plot.figure_id);hold on;
    [~,~]=pcontour(lp3.V,double(CC),plot.domain,'c');
    solU = []; optimal_vc = []; lp1iter = 0;
    %%
    while abs(double(CC)-double(C0)) >= 1e-6
        lp1iter = lp1iter + 1;
        if lp1iter ~= 0
            C0 = CC;
        end
        [solL,lp4.feas]= sos_optimal_v2(sym.f,sym.gg,lp4.u,lp4.au,lp4.V,C0);
        if lp4.feas == 0
            break
        elseif isempty(p2s(solL))
            fprintf('196 The auxilirary of optimal LF does not exist.======\n');
        else
            fprintf('198 The sublevel set of optimal LF move on.======\n');
        end
        [CC,lp4.feas,solu1,solu2] = sos_optimal_v3(sym.f,sym.gg,lp4.u,lp4.au,lp4.V,sym.us,dom,solL,ccc,plot.figure_id);
        optimal_vc = [optimal_vc; double(CC)]
        solU = [solU;[solu1,solu2]];
        if lp4.feas == 0
            break
        end
    end
    figure(plot.figure_id);hold on;
    [~,~]=pcontour(lp4.V,optimal_vc(end),plot.domain,'k');
    %% Start to compute the control barrier function
    c_b = optimal_vc(end);
    sol_B = c_b - N_Lya(end);
    %% Hyperparameters of the SOSP @ CBF LP5
    %     lp5.u = 4;
    %     lp5.h = 4;
    %     lp5.us = 8;
    %     lp5.au = 8;
    %     lp5.gamma = 0;
    %%
    lp5.V = N_Lya(end);
    solh = sol_B;
    trace_Q1 = 1; trace_Q = 0;
    lp5.feas = 1; j = 0;
    %%
    figure(plot.figure_id+1);clf;hold on;
    [~,~]=pcontour(sym.us(1),0,plot.domain,'r');
    [~,~]=pcontour(sym.us(2),0,plot.domain,'r');
    if length(sym.us) == 3
        [~,~]=pcontour(sym.us(3),0,plot.domain,'r');
    elseif length(sym.us) == 4
        [~,~]=pcontour(sym.us(3),0,plot.domain,'r');
        [~,~]=pcontour(sym.us(4),0,plot.domain,'r');
    else
        fprintf('233 The constraint number does not match.======\n');
    end
    record_Q = [1];
    trace_Q = 1.1;
    while abs(double(trace_Q)-double(record_Q(end)))>=1e-3
        j = j+1
        record_Q = trace_Q
        [SOLu1,SOLu2,SOL1,SOL2,lp5.feas] = sos_function_1(sym.f,lp5.u,lp5.au,solh,lp5.V,lp5.gamma,sym.gg);
        if lp5.feas == 0
            fprintf('242 Optimal CBF controller [End] can not find.======\n');
            break
        end
        Control_lp5 = [Control_lp5; [SOLu1 SOLu2]];
        [solh,trace_Q,lp5.feas]=sos_function_2(j,sym.f,b_k_h,SOLu1,SOLu2,SOL1,SOL2,lp5.gamma,lp5.V,sym.us,dom,sym.gg,lp5.us,plot.figure_id+3);
        TRACE_lp5 = [TRACE_lp5; double(trace_Q)];
        Barrier_lp5 = [Barrier_lp5; solh];
        if lp5.feas == 0
            fprintf('250 Optimal Control Barrier function can not find.======\n');
            break
        end
    end
    %%
    figure(plot.figure_id);hold on;
    if ~isempty(Barrier_lp5)
        [~,~]=pcontour(Barrier_lp5(end),0,plot.domain,'b');
    end
    k_iter = 5;
end
end