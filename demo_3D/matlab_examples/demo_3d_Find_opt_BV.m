clear;tic;
pvar x1 x2 x3 u htol epsi;
x = [x1;x2;x3];
% load test_1125_3_3D_opt_Lya;
% load test_1126_3_3D_opt_Lya;
% load test_1126_3_3D_to_OPT_Lya;
% load test_1126_4_3D_to_OPT_Lya;
load test_1127_opt_Lya_4424_444_44_result_output;
dom = 10; domain = [-dom dom -dom dom -dom dom];
figure_id = 311;
%%
C1 = (x1-3)^2+(x2-0)^2+(x3-0)^2-2;
C2 = (x1+3)^2+(x2+3)^2+(x3-3)^2-3;
C3 = (x1-0)^2+(x2-3)^2+(x3+0)^2-3;
C4 = (x1-3)^2+(x2-0)^2+(x3+3)^2-3;
C = [C2;C3;C4];
V00 = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
CC0 = 2.584683714740699;
%%
Compute = {};
%%
% for i = [12;14;15;18;22;24;25;26;27;28;29;30;31]
% for i = 1:length(A)
for i = 15
    %%
    f = [-0.16211179709695037165964324641892*x1^4+0.49031898059485488031675707268867*x1^3-0.80741059860165544525334054266171*x1^2-0.67918469453797989758096302163419*x1-0.000031796033448678815965058458425929*x2^4+0.00049811603934539429219124917480599*x2^3-0.017884574789259536503616132563366*x2^2+0.059003602596920091960530641017613*x2+0.14736298331076760903535216584714*x3^4-0.22234685217253638556123007674614*x3^3+0.14008297734444075111071015271591*x3^2+0.34076230832989312657943514750514*x3-1.4865840251194628709408125437986
        -x2*x1^3-x2
        1.876321954439673866943394386908*x1^4-1.8803947547684580765547934788628*x1^3-x1^2*x3-0.3521313244010295662178577913437*x1^2-0.58570314276461321600919518459705*x1+0.0008606977977094753011130801034767*x2^4-0.0047246616639717428295930368165045*x2^3+0.0073970075005580443808228530144788*x2^2-0.39277632595769917944750204696902*x2-0.63399852647351906398398568853736*x3^4-2.3526877462367330462456038731034*x3^3-0.18912599086086179234200699283974*x3^2-1.2624646738177363047839207865763*x3+1.2361363889678920191528277428006
        ];
    gg = [1;1;1];
    %%
    trace_Q1 = 1; trace_Q = 0;
    mm = 0; kk = 1;
    %%
    %     k_u = 4; k_h = 4; L_us = 6; L_au = 4; gamma = 0;
    k_u = 4; k_h = 4; L_us = 2; L_au = 2; gamma = 0;
    %%
    %     for j = 1:14
    %         if j == 1
    %             % Get
    %             k_u = 4; k_h = 4; L_us = 4; L_au = 4; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 2
    %             k_u = 4; k_h = 4; L_us = 6; L_au = 8; gamma = 0;
    %             figure_id = figure_id + 10;
    %         elseif j == 3
    %             k_u = 4; k_h = 4; L_us = 8; L_au = 6; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 4
    %             k_u = 4; k_h = 4; L_us = 8; L_au = 8; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 5
    %             k_u = 4; k_h = 4; L_us = 6; L_au = 6; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 6
    %             k_u = 4; k_h = 4; L_us = 4; L_au = 6; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 7
    %             % Get
    %             k_u = 4; k_h = 4; L_us = 6; L_au = 4; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 8
    %             k_u = 6; k_h =4; L_us = 6; L_au = 4; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 9
    %             k_u = 6; k_h = 6; L_us = 6; L_au = 6; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 10
    %             k_u = 6; k_h = 6; L_us = 8; L_au = 6; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 11
    %             k_u = 6; k_h = 6; L_us = 6; L_au = 8; gamma = 0;
    %             figure_id = figure_id + 100;
    %         elseif j == 12
    %             k_u = 6; k_h = 6; L_us = 8; L_au = 8; gamma = 0;
    %             figure_id = figure_id + 100;
    %         else
    %             k_u = 8; k_h = 8; L_us = 8; L_au = 8; gamma = 0;
    %             figure_id = figure_id + 100;
    %         end
    %%
    figure(figure_id);clf;hold on;
    %%
    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);
    B_controller = A(i,4);
    V = A_V(i);
    C0 = A_C0(i);
    sol_B = C0 - V; solh = sol_B;
    %     phB= patch(pcontour3(B_controller,0,domain,'B')); set(phB,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',1);
    %%
    TRACE = []; Barrier = []; Control = []; iter = 0;
    %%
    while 1
        %%
        %         B1 = -9.161325193733466*x1^4-0.03439578029383424*x1^3*x2+0.05334550884420283*x1^3*x3+16.284547859344*x1^2*x2^2+0.0635207097935428*x1^2*x2*x3+1.671630410877847*x1^2*x3^2-0.03374348314766387*x1*x2^3+0.06596074193275561*x1*x2^2*x3+0.03421307872205565*x1*x2*x3^2-0.03777639021097199*x1*x3^3-9.178911093916067*x2^4+0.05280503556005673*x2^3*x3+1.661285535205469*x2^2*x3^2-0.04117062123967289*x2*x3^3-1.866306920312326*x3^4-0.00814177614973843*x1^3-0.01242219063700423*x1^2*x2-0.008978724863385771*x1^2*x3-0.02639232438017272*x1*x2^2-0.00530455804716668*x1*x2*x3+0.02386596853139026*x1*x3^2-0.02161135252213946*x2^3+0.00254072715075745*x2^2*x3+0.03869785494639697*x2*x3^2+0.002202029248068501*x3^3-0.07773539805278774*x1^2+0.009715359497566189*x1*x2-0.04940480294142099*x1*x3-0.06280416965745615*x2^2-0.05560402831825616*x2*x3-0.3686872166926973*x3^2-1.356091452933523e-05*x1+2.817871385497815e-06*x2+4.698187562390626e-05*x3+14.06958796419195;
        B1 = A(i,4);
        inV = patch(pcontour3(B1,0,domain,'b'));              % Plot the original Lyapunov sublevel set
        set(inV,'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','-','LineWidth',0.6); hold on;
        
        ph2= patch(pcontour3(C2,0,domain,'c')); set(ph2, 'FaceColor', 'none', 'EdgeColor', 'k' );
        ph3= patch(pcontour3(C3,0,domain,'c')); set(ph3, 'FaceColor', 'none', 'EdgeColor', 'k' );
        ph4= patch(pcontour3(C4,0,domain,'c')); set(ph4, 'FaceColor', 'none', 'EdgeColor', 'k' );
        phV0= patch(pcontour3(V,double(C0),domain,'G')); set(phV0,'EdgeAlpha',1, 'FaceColor', 'none', 'EdgeColor', 'g','LineStyle','-','LineWidth',1);
        %%
        iter = iter+1
        record_Q = trace_Q
        %%
        [SOLu1,SOLu2,SOLu3,SOL1,SOL2,kk] = sos_function_1_3D(f,k_u,L_au,solh,V,gamma,gg);
        if kk == 0
            break
        end
        Control = [Control; [SOLu1 SOLu2 SOLu3]];
        %%
        [solh,trace_Q,kk] = sos_function_2_3D(iter,f,k_h,SOLu1,SOLu2,SOLu3,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id);
        if kk == 0
            break
        end
        TRACE = [TRACE; double(trace_Q)];
        Barrier = [Barrier; solh];
    end
    view(-121,-10);
    figure_id = figure_id+1;
    figure(figure_id);clf;
    Compute{i} = {Barrier};
    %     end
end