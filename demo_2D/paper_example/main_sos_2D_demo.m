clear;
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
% f = [x2; -x1-x2*(1-x1^2)];
f = [x2-x1;
    0.46382156330341729589940137480476*x1^4-0.21219821363526560116570580989732*x1^3+0.35000677410250802565647474138031*x1^2+x1*x2-0.27258976376810471290805063896793*x1+0.15407267311901634565529661813343*x2^4+0.94268273772394006737584959410015*x2^3+4.1277331008261466394060335005634*x2^2-1.7098389065959235244562819389103*x2+0.011510115809394316777058975276304
    ];
gg = [1; 1];
sys = [f(1)+gg(1)*u1; f(2)+gg(2)*u2];
%%
% C1 = (x1+3)^2+(x2-4)^2-1;
% C2 = (x1+3)^2+(x2+4)^2-1;
% C3 = (x1-3)^2+(x2-2)^2-1;
C1 = (x1+8)^2+(x2-0)^2-4;
C2 = (x1-8)^2+(x2+0)^2-4;
C3 = (x1-0)^2+(x2-8)^2-4;
C4 = (x1-0)^2+(x2+8)^2-4;
C = [C1;C2;C3;C4];
% C = [C1;C2;C3];
% V0 = x1^2+x1*x2+x2^2;
% V0 = 4*x1^4+2*x2^4+2*x1^2*x2^2+4*x1^2+2*x2^2+1*x1*x2;
V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; 
c0 = 1; cc = 1.1; epsi = 1e-6; epsi2 = 0.196;
figure_id = 111;
%%
dom = 15; domain = [-dom dom -dom dom];
figure(figure_id);clf;hold on;
[~,~]=pcontour(C(1),0,domain,'r');
[~,~]=pcontour(C(2),0,domain,'r');
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'r');
elseif length(C) == 4
    [~,~]=pcontour(C(3),0,domain,'r');
    [~,~]=pcontour(C(4),0,domain,'r');
else
    fprintf('The constraint number does not match.======\n');
end
[~,~]=pcontour(V0,c0,domain,'g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hyperparameter of the SOSP @ CLF
solU = []; v_c = []; iter = 1; kk = 1;
% v_k_u = 4; 
v_k_u = 2;
v_k_l = 4;
%% Start to compute the sublevel set of a Lyapunov function
% while abs(double(cc)-double(c0)) >= epsi
while double(cc)-double(c0) >= epsi
    iter = iter + 1;
    if iter ~= 1
        c0 = cc;
    end
    [solu1,solu2,solL,kk]= sos_function_v(f,gg,v_k_u,v_k_l,V0,c0);
    if kk == 0
        break
    end
    [cc,kk,solu1,solu2] = sos_function_v2(f,gg,v_k_u,v_k_l,V0,C,dom,solL,figure_id);
    double(cc)
    v_c = [v_c; double(cc)];
    solU = [solU;[solu1,solu2]];
    if kk == 0
        break
    end
end
hold on;
[~,~]=pcontour(V0,max(double(v_c)),domain,'b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start to compute the control barrier function
c_b = max(double(v_c));
sol_B = c_b - V0;
solh = sol_B;
%% Hyperparameters of the SOSP @ CBF
% b_k_u = 6; b_k_h = 4;
% L_us = 8; L_au = 8;
% gamma = 0;
b_k_u = 4; b_k_h = 4;
L_us = 8; L_au = 8;
gamma = 0;
kk = 1; j = 0;
TRACE = [];
Barrier = [];
Control = [];
%%
figure(figure_id+1);clf;hold on;
[~,~]=pcontour(C(1),0,domain,'r');
[~,~]=pcontour(C(2),0,domain,'r');
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'r');
elseif length(C) == 4
    [~,~]=pcontour(C(3),0,domain,'r');
    [~,~]=pcontour(C(4),0,domain,'r');
else
    fprintf('The constraint number does not match.======\n');
end
record_Q = [1000];
trace_Q = 10001;
while abs(double(trace_Q)-double(record_Q(end)))>=epsi2
    j = j+1
    record_Q = [record_Q; trace_Q];
    record_Q
    [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1(f,b_k_u,L_au,solh,V0,gamma,gg);
    if kk == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2]];
    [solh,trace_Q,kk]=sos_function_2(j,f,b_k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V0,C,dom,gg,L_us,figure_id+1);
    if kk == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
end
%%
figure(figure_id);hold on;
if ~isempty(Barrier)
    [~,~]=pcontour(Barrier(end),0,domain,'m');
else
    fprintf('Barrier function can not find.======\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1
    %% Hyperparameters of the SOSP @ CBF -> V
    V_us = 4; V_au = 4; V_degree = 4;
    gamma = 0;
    kk = 1; OO = 0;
    %%
    i_count = min(length(Control),length(Barrier));
    if ~isempty(i_count)
        u1 = Control(i_count,1);
        u2 = Control(i_count,2);
        B = Barrier(i_count);
    else
        u1 = Control(i_count,1);
        u2 = Control(i_count,2);
        B = 0;
        fprintf('Barrier function can not find.======\n');
    end

    dom_2 = 100; domain_2 = [-dom_2 dom_2 -dom_2 dom_2];
    N_Lya = [];
    %%
    [V, kk] = sos_optimal_V1(f,gg,B,u1,u2,V_au,V_us,V_degree,C,gamma);
    %%
    if kk == 0
        fprintf('Suitable Lyapunov function can not find.======\n');
    end
    N_Lya = [N_Lya;V];
    figure(figure_id+2);clf;hold on;
    cc = [1,2,3];
    for i = 0:16
        k_u_V = ['r','g','b','m','c','k','y'];
        if mod(i,7) == 0
            [~,~]=pcontour(V,cc(1)*i,domain_2,k_u_V(7)); hold on;
            [~,~]=pcontour(V,cc(2)*i,domain_2,k_u_V(7)); hold on;
            [~,~]=pcontour(V,cc(3)*i,domain_2,k_u_V(7)); hold on;
        else
            [~,~]=pcontour(V,cc(1)*i,domain_2,k_u_V(mod(i,7))); hold on;
            [~,~]=pcontour(V,cc(2)*i,domain_2,k_u_V(mod(i,7))); hold on;
            [~,~]=pcontour(V,cc(3)*i,domain_2,k_u_V(mod(i,7))); hold on;
        end
    end
    Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
    for i = 0:16
        k_u_V = ['r','g','b','m','c','k','y'];
        if mod(i,7) == 0
            [~,~]=pcontour(Vdot,0,domain_2,k_u_V(7)); hold on;             % Plot the original Lyapunov sublevel set
        else
            [~,~]=pcontour(Vdot,0,domain_2,k_u_V(mod(i,7))); hold on;             % Plot the original Lyapunov sublevel set
        end
    end
    %% TEST FOR Sublevel Set
    [a1,b1] = coeffs(p2s(V));
    ccc = double(vpa(a1(end)));
    C0 = ccc; cc = ccc+0.1;
    k_u_V = 4; k_l_au = 4; kk = 1;
    %%
    figure(figure_id);hold on;
    [~,~]=pcontour(V,double(cc),domain,'c');
    solU = []; v_c = []; iter = 0;
    %%
    while abs(double(cc)-double(C0)) >= epsi
        iter = iter + 1;
        if iter ~= 0
            C0 = cc;
        end
        [solu1,solu2,solL,kk]= sos_optimal_v2(f,gg,k_u_V,k_l_au,V,C0);
        if kk == 0
            break
        end
        [cc,kk,solu1,solu2] = sos_optimal_v3(f,gg,k_u_V,k_l_au,V,C,dom,solL,ccc,figure_id);
        v_c = [v_c; double(cc)]
        solU = [solU;[solu1,solu2]];
        if kk == 0
            figure(figure_id);hold on;
            [~,~]=pcontour(V,v_c(end),domain,'b');
            break
        end
    end
    
    %% Start to compute the control barrier function
    c_b = v_c(end);
    sol_B = c_b - N_Lya(end);
    figure(figure_id);hold on;
    [~,~]=pcontour(sol_B,0,domain,'k');
    V = N_Lya(end);
    %% Hyperparameters of the SOSP @ CBF
    solh = sol_B;
    b_k_u = 4; b_k_h = 4;
    L_us = 8; L_au = 8;
    gamma = 1;
    trace_Q1 = 1; trace_Q = 0;
    kk = 1; j = 0;
    TRACE = [];
    Barrier = [];
    Control = [];
    %%
    figure(figure_id+1);clf;hold on;
    [~,~]=pcontour(C(1),0,domain,'r');
    [~,~]=pcontour(C(2),0,domain,'r');
    if length(C) == 3
        [~,~]=pcontour(C(3),0,domain,'r');
    elseif length(C) == 4
        [~,~]=pcontour(C(3),0,domain,'r');
        [~,~]=pcontour(C(4),0,domain,'r');
    else
        fprintf('The constraint number does not match.======\n');
    end
    record_Q = [0.0001];
    trace_Q = 0.00011;
    while abs(double(trace_Q)-double(record_Q(end)))>=5e-1
        j = j+1
        record_Q = trace_Q
        [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1(f,b_k_u,L_au,solh,V,gamma,gg);
        if kk == 0
            break
        end
        Control = [Control; [SOLu1 SOLu2]];
        [solh,trace_Q,kk]=sos_function_2(j,f,b_k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id+1);
        if kk == 0
            fprintf('Optimal Control Barrier function can not find.======\n');
            break
        end
        TRACE = [TRACE; double(trace_Q)];
        Barrier = [Barrier; solh];
    end
    %%
    figure(figure_id);hold on;
    [~,~]=pcontour(Barrier(end),0,domain,'k');
end