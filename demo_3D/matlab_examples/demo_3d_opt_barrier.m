clear;tic;
pvar x1 x2 x3 u htol epsi;
x = [x1;x2;x3];
%%
f = [-0.16211179709695037165964324641892*x1^4+0.49031898059485488031675707268867*x1^3-0.80741059860165544525334054266171*x1^2-0.67918469453797989758096302163419*x1-0.000031796033448678815965058458425929*x2^4+0.00049811603934539429219124917480599*x2^3-0.017884574789259536503616132563366*x2^2+0.059003602596920091960530641017613*x2+0.14736298331076760903535216584714*x3^4-0.22234685217253638556123007674614*x3^3+0.14008297734444075111071015271591*x3^2+0.34076230832989312657943514750514*x3-1.4865840251194628709408125437986
    -x2*x1^3-x2
    1.876321954439673866943394386908*x1^4-1.8803947547684580765547934788628*x1^3-x1^2*x3-0.3521313244010295662178577913437*x1^2-0.58570314276461321600919518459705*x1+0.0008606977977094753011130801034767*x2^4-0.0047246616639717428295930368165045*x2^3+0.0073970075005580443808228530144788*x2^2-0.39277632595769917944750204696902*x2-0.63399852647351906398398568853736*x3^4-2.3526877462367330462456038731034*x3^3-0.18912599086086179234200699283974*x3^2-1.2624646738177363047839207865763*x3+1.2361363889678920191528277428006
    ];
gg = [1;1;1];
% input = [gg(1)*u1;gg(2)*u2;gg(2)*u3];
V = 10*x1^4+1*x2^4+20*x3^4+2*x1^2*x2^2-4*x3^2*x2^2+3*x1^2*x3^2;
C0 = 96.811595465770495;
% input = [gg(1)*u1;gg(2)*u2;gg(3)*u3];
%%
C1 = (x1+4)^2+(x2-6)^2+(x3+2)^2-4;
C2 = (x1+3)^2+(x2+4)^2+(x3+4)^2-4;
C3 = (x1-4)^2+(x2-0)^2+(x3-0)^2-5;
C4 = (x1+4)^2+(x2-2)^2+(x3-4)^2-5;
C = [C1;C2;C3;C4];
trace_Q1 = 1; trace_Q = 0;
mm = 0; kk = 1; i = 0;
%%
sol_B = C0 - V;
solh = sol_B;
%%
k_u = 4;
k_h = 4;
L_us = 4;
L_au = 4;
gamma = 0;
%%
figure_id = 22;
figure(figure_id);clf;hold on;
dom = 10; domain = [-dom dom -dom dom -dom dom];
phV0= patch(pcontour3(V,double(C0),domain,'c')); set(phV0,'EdgeAlpha',0.1, 'FaceColor', 'none', 'EdgeColor', 'b' );
ph1= patch(pcontour3(C1,0,domain,'c')); set(ph1, 'FaceColor', 'none', 'EdgeColor', 'k' );
ph2= patch(pcontour3(C2,0,domain,'c')); set(ph2, 'FaceColor', 'none', 'EdgeColor', 'k' );
ph3= patch(pcontour3(C3,0,domain,'c')); set(ph3, 'FaceColor', 'none', 'EdgeColor', 'k' );
ph4= patch(pcontour3(C4,0,domain,'c')); set(ph4, 'FaceColor', 'none', 'EdgeColor', 'k' );
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(-30,20);

TRACE = []; Barrier = []; Control = [];
%%
while 1
    i = i+1
    record_Q = trace_Q
    %%
    [SOLu1,SOLu2,SOLu3,SOL1,SOL2,kk] = sos_function_1_3D(f,k_u,L_au,solh,V,gamma,gg);
    if kk == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2 SOLu3]]
    %%
    [solh,trace_Q,kk] = sos_function_2_3D(i,f,k_h,SOLu1,SOLu2,SOLu3,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id);
    if kk == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
end