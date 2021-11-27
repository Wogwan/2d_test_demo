clear;
% load test_1126_2_3D;
load test_1127_opt_barrier_4424_result_output;
pvar x1 x2 x3 u1 u2 u3 htol epsi;
format long
x = [x1;x2;x3];
%%
f = [-0.16211179709695037165964324641892*x1^4+0.49031898059485488031675707268867*x1^3-0.80741059860165544525334054266171*x1^2-0.67918469453797989758096302163419*x1-0.000031796033448678815965058458425929*x2^4+0.00049811603934539429219124917480599*x2^3-0.017884574789259536503616132563366*x2^2+0.059003602596920091960530641017613*x2+0.14736298331076760903535216584714*x3^4-0.22234685217253638556123007674614*x3^3+0.14008297734444075111071015271591*x3^2+0.34076230832989312657943514750514*x3-1.4865840251194628709408125437986
    -x2*x1^3-x2
    1.876321954439673866943394386908*x1^4-1.8803947547684580765547934788628*x1^3-x1^2*x3-0.3521313244010295662178577913437*x1^2-0.58570314276461321600919518459705*x1+0.0008606977977094753011130801034767*x2^4-0.0047246616639717428295930368165045*x2^3+0.0073970075005580443808228530144788*x2^2-0.39277632595769917944750204696902*x2-0.63399852647351906398398568853736*x3^4-2.3526877462367330462456038731034*x3^3-0.18912599086086179234200699283974*x3^2-1.2624646738177363047839207865763*x3+1.2361363889678920191528277428006
    ];
gg = [1;1;1];
sys = [f(1)+gg(1)*u1; f(2)+gg(2)*u2; f(3)+gg(3)*u3];
%%
dom = 10; domain = [-dom dom -dom dom -dom dom];
C1 = (x1-3)^2+(x2-0)^2+(x3-0)^2-2;
C2 = (x1+3)^2+(x2+3)^2+(x3-3)^2-3;
C3 = (x1-0)^2+(x2-3)^2+(x3+0)^2-3;
C4 = (x1-3)^2+(x2-0)^2+(x3+3)^2-3;
C = [C2;C3;C4];
V0 = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
CC0 = 2.584683714740699;
figure_id = 211;
A_V = [];
A_C0 = [];
for iter = 1:length(A)
    %% Start to compute the control barrier function
    u1 = A(iter,1);
    u2 = A(iter,2);
    u3 = A(iter,3);
    B = A(iter,4);
    %%
    c0 = 1; cc = 1.1; epsi = 1e-6; N_Lya = [];
    %%
    figure(figure_id);clf;hold on;
    %     ph1= patch(pcontour3(C1,0,domain,'c')); set(ph1, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph2= patch(pcontour3(C2,0,domain,'c')); set(ph2, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph3= patch(pcontour3(C3,0,domain,'c')); set(ph3, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph4= patch(pcontour3(C4,0,domain,'c')); set(ph4, 'FaceColor', 'none', 'EdgeColor', 'k' );
    phV0= patch(pcontour3(V0,double(CC0),domain,'G')); set(phV0,'EdgeAlpha',1, 'FaceColor', 'none', 'EdgeColor', 'g' );
    phB= patch(pcontour3(B,0,domain,'B')); set(phB,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',1);
    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(-100,0);
    %% Hyperparameters of the SOSP @ CBF -> V
    V_us = 4; V_au = 4; V_degree = 4; gamma = 0; 
    k_u_V = 4; k_l_au = 4;
    % V_us = 6; V_au = 6; V_degree = 6; gamma = 0; k_u_V = 6; k_l_au = 6;
    % V_us = 8; V_au = 8; V_degree = 8; gamma = 0; k_u_V = 8; k_l_au = 8;
    %%
    kk = 1; OO = 0;
    %%
    % [V, kk] = sos_optimal_V1_3D(f,gg,B,u1,u2,u3,V_au,V_us,V_degree,C,gamma);
    [V, kk] = sos_optimal_V1_3D_V(f,gg,B,u1,u2,u3,V_au,V_us,V_degree,C,gamma);
    %%
    if kk == 0
        fprintf('Suitable Lyapunov function can not find.======\n');
    end
    N_Lya = [N_Lya;V];
    %% TEST FOR Sublevel Set
    [a1,b1] = coeffs(p2s(V));
    ccc = double(vpa(a1(end)));
    C0 = ccc; cc = ccc+1;solU = []; v_c = []; iter = 0;
    %%
    figure(figure_id);hold on;
    %     phV1= patch(pcontour3(V,double(cc),domain,'G')); set(phV1,'EdgeAlpha',1, 'FaceColor', 'none', 'EdgeColor', 'r' );
    %%
    while double(cc)-double(C0) >= epsi
        iter = iter + 1;
        if iter ~= 0
            C0 = cc;
        end
        [solL,kk]= sos_optimal_v2_3D(f,gg,k_u_V,k_l_au,V,cc);
        if kk == 0
            break
        end
        [cc,kk,solu1,solu2,solu3] = sos_optimal_v3_3D(f,gg,k_u_V,k_l_au,V,C,dom,solL,ccc,figure_id);
        v_c = [v_c; double(cc)]
        solU = [solU;[solu1,solu2,solu3]];
    end
    %% Start to compute the control barrier function
    c_b = max(double(v_c));
    V = N_Lya(end);
    A_V = [A_V;V];
    A_C0 = [A_C0; double(c_b)];
    sol_B = c_b - N_Lya(end);
    hold on;
    phsol_B= patch(pcontour3(sol_B,0,domain,'B')); set(phsol_B,'EdgeAlpha',0.7,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','-','LineWidth',1);
    
    figure_id = figure_id + 1;
    figure(figure_id);clf;
end