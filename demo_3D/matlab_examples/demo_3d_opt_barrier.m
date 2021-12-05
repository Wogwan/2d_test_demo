clear;tic;
pvar x1 x2 x3 u htol epsi;
x = [x1;x2;x3];
%%
f = [-0.16211179709864792508611230914539*x1^4+0.46967680241313530098423711933719*x1^3-0.80741059859268821864900068453951*x1^2-0.52420770456533548609101558213297*x1-0.000031796033448678815965058458425929*x2^4+0.00049811603934539429219124917480599*x2^3-0.017884574789259536503616132563366*x2^2+0.059003602596920091960530641017613*x2+0.14736298331076760903535216584714*x3^4-0.22234685217253638556123007674614*x3^3+0.14008297734444075111071015271591*x3^2+0.34076230832989312657943514750514*x3-1.4865840251312778030530597727999
    -x2*x1^3-x2
    1.876321954439673866943394386908*x1^4-1.8803947547684580765547934788628*x1^3-x1^2*x3-0.3521313244010295662178577913437*x1^2-0.58570314276461321600919518459705*x1+0.0008606977977094753011130801034767*x2^4-0.0047246616639717428295930368165045*x2^3+0.0073970075005580443808228530144788*x2^2-0.39277632595769917944750204696902*x2-6.3527586089398613289347395038931*x3^4+9.1026400980830752818206974552595*x3^3+5.2275118144867558367394622109714*x3^2-10.044858323825100132609122738359*x3+0.61281374087452245014162599545671
    ];
gg = [1;1;1];
%%
C2 = (x1+4)^2+(x2+4)^2+(x3-4)^2-4;
C3 = (x1-0)^2+(x2-4)^2+(x3+0)^2-4;
C4 = (x1-4)^2+(x2-0)^2+(x3+4)^2-6;
C = [C2;C3;C4];
trace_Q1 = 1; trace_Q = 0;
mm = 0; kk = 1;
%%
% V = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2; % C0 = 16.00000001852458;
V = x1^4+x1^3*x2+x1^2*x2^2+x1*x2^3+x2^4+x1^3*x3+x1^2*x2*x3+x1*x2^2*x3+x2^3*x3+x1^2*x3^2+x1*x2*x3^2+x2^2*x3^2+x1*x3^3+x2*x3^3+x3^4+x1^3+x1^2*x2+x1*x2^2+x2^3+x1^2*x3+x1*x2*x3+x2^2*x3+x1*x3^2+x2*x3^2+x3^3+x3^2+x1^2+x1*x2+x2^2+x1*x3+x2*x3+x1+x2+x3+1; 
C0 = 24.92234264283737;
%%
sol_B = C0 - V; solh = sol_B;
%%
k_u = 4; k_h = 4; L_us = 4; L_au = 4; gamma = 0;
%%
figure_id = 1;
figure(figure_id);clf;hold on;
dom = 10; domain = [-dom dom -dom dom -dom dom];
%%
TRACE = []; Barrier = []; Control = []; i = 1;
%%
% for i = 1:10
while 1
    record_Q = trace_Q
    phV0= patch(pcontour3(V,double(C0),domain,'c')); set(phV0,'EdgeAlpha',0.1, 'FaceColor', 'none', 'EdgeColor', 'b' );
    ph2= patch(pcontour3(C2,0,domain,'c')); set(ph2, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph3= patch(pcontour3(C3,0,domain,'c')); set(ph3, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph4= patch(pcontour3(C4,0,domain,'c')); set(ph4, 'FaceColor', 'none', 'EdgeColor', 'k' );
    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(45,8);
    %%
    [SOLu1,SOLu2,SOLu3,SOL1,SOL2,kk] = sos_function_1_3D(f,k_u,L_au,solh,V,gamma,gg);
    if kk == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2 SOLu3]];
    %%
    [solh,trace_Q,kk] = sos_function_2_3D(i,f,k_h,SOLu1,SOLu2,SOLu3,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id);
    if kk == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
    if i ~= 15
        figure_id = figure_id+1;
        figure(figure_id);clf;
    end
end

A = [];
for iter = 1:length(Barrier)
    A = [A; [Control(iter,:) Barrier(iter)]];
end
%%
fprintf('Permissive B(x) is \n%s \n\n',char(vpa(p2s(A(end,4)))));
fprintf('Control Input u1(x) is \n%s \n\n',char(vpa(p2s(A(end,1)))));
fprintf('Control Input u2(x) is \n%s \n\n',char(vpa(p2s(A(end,2)))));
fprintf('Control Input u3(x) is \n%s \n\n',char(vpa(p2s(A(end,3)))));