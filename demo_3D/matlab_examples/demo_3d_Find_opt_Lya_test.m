clear; load BC_test_new_Lya;
pvar x1 x2 x3 u1 u2 u3 htol epsi;
format long
x = [x1;x2;x3];
%%
V_Sel = [];
C0_Sel = [];
for i = 3:length(A)
    f = [-0.16211179709695037165964324641892*x1^4+0.49031898059485488031675707268867*x1^3-0.80741059860165544525334054266171*x1^2-0.67918469453797989758096302163419*x1-0.000031796033448678815965058458425929*x2^4+0.00049811603934539429219124917480599*x2^3-0.017884574789259536503616132563366*x2^2+0.059003602596920091960530641017613*x2+0.14736298331076760903535216584714*x3^4-0.22234685217253638556123007674614*x3^3+0.14008297734444075111071015271591*x3^2+0.34076230832989312657943514750514*x3-1.4865840251194628709408125437986
        -x2*x1^3-x2
        1.876321954439673866943394386908*x1^4-1.8803947547684580765547934788628*x1^3-x1^2*x3-0.3521313244010295662178577913437*x1^2-0.58570314276461321600919518459705*x1+0.0008606977977094753011130801034767*x2^4-0.0047246616639717428295930368165045*x2^3+0.0073970075005580443808228530144788*x2^2-0.39277632595769917944750204696902*x2-0.63399852647351906398398568853736*x3^4-2.3526877462367330462456038731034*x3^3-0.18912599086086179234200699283974*x3^2-1.2624646738177363047839207865763*x3+1.2361363889678920191528277428006
        ];
    gg = [1;1;1];
    sys = [f(1)+gg(1)*u1; f(2)+gg(2)*u2; f(3)+gg(3)*u3];
    c0 = 1; cc = 1.1; epsi = 1e-6; figure_id = 211; N_Lya = [];
    %%
    C1 = (x1+4)^2+(x2-6)^2+(x3+2)^2-4;
    C2 = (x1+3)^2+(x2+4)^2+(x3+4)^2-4;
    C3 = (x1-4)^2+(x2-0)^2+(x3-0)^2-5;
    C4 = (x1+4)^2+(x2-2)^2+(x3-4)^2-5;
    %     C = [C1;C2;C3;C4];
    C = [C2;C3;C4];
    %%
    V0 = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
    C0 = 9.681159889041236;
    %% Start to compute the control barrier function 24
    u1 = A(i,1);
    u2 = A(i,2);
    u3 = A(i,3);
    B = A(i,4);
    %%
    figure(figure_id);clf;hold on;
    dom = 10; domain = [-dom dom -dom dom -dom dom];
    %     ph1= patch(pcontour3(C1,0,domain,'c')); set(ph1, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph2= patch(pcontour3(C2,0,domain,'c')); set(ph2, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph3= patch(pcontour3(C3,0,domain,'c')); set(ph3, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph4= patch(pcontour3(C4,0,domain,'c')); set(ph4, 'FaceColor', 'none', 'EdgeColor', 'k' );
    phV0= patch(pcontour3(V0,double(C0),domain,'G')); set(phV0,'EdgeAlpha',1, 'FaceColor', 'none', 'EdgeColor', 'g' );
    phB= patch(pcontour3(B,0,domain,'B')); set(phB,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',1);
    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(-30,20);
    %%
    figure(figure_id+100);clf;hold on;
    phV0= patch(pcontour3(V0,double(C0),domain,'G')); set(phV0,'EdgeAlpha',1, 'FaceColor', 'none', 'EdgeColor', 'g' );
    %% Hyperparameters of the SOSP @ CBF -> V
    V_us = 4; V_au = 6; V_degree = 4; k_u_V = 4; k_l_au = 4; gamma = 0;
    % V_us = 6; V_au = 6; V_degree = 6; gamma = 0; k_u_V = 6; k_l_au = 6;
    % V_us = 8; V_au = 8; V_degree = 8; gamma = 0; k_u_V = 8; k_l_au = 8;
    %%
    kk = 1; OO = 0;
    %%
    [V, kk] = sos_optimal_V1_3D(f,gg,B,u1,u2,u3,V_au,V_us,V_degree,C,gamma);
    %     [V, kk] = sos_optimal_V1_3D_Update(f,gg,B,u1,u2,u3,V_au,V_us,V_degree,C,gamma);
    %% Start to compute the control barrier function
    phsol_B= patch(pcontour3(V,0,domain,'B')); 
    set(phsol_B,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'c','LineStyle','-','LineWidth',1);
    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(-30,20);
    figure(figure_id);hold on;
    phsol_B= patch(pcontour3(V,0,domain,'B')); 
    set(phsol_B,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','-','LineWidth',1);
    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(-30,20);
    figure_id = figure_id + 1;

end

B = [];
for i = 1:length(V_Sel)
    B = [B; [V_Sel(i) C0_Sel(i)]];
end