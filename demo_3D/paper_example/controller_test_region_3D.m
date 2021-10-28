clear;
pvar x1 x2 x3 u1 u2 htol epsi;
x = [x1;x2;x3];
%%%%%%%%%%%%%%%%%%%%%%%
%%
% f = [x2-x3^2; x3-x1^2; -x1-2*x2-x3+x2^3];
f = [x3^2+x2
    0.20168850953378745366473268063601*x1^4+0.063529983054143205658580619337485*x1^3-1.5117281788498260486353075293664*x1^2+0.11050308587243914680550528927443*x1-0.010789899544848068069224922282956*x2^4-0.0023351458555991655777206439381644*x2^3-0.050856189012352337464051288407063*x2^2-0.011697431254263418801131457769316*x2-0.0050275698408848868564691159122049*x3^4-0.072414399333034673578168849417125*x3^3+0.043478883248082855761396103844163*x3^2+1.1021325793385576968796968344577*x3-0.41904249821184156449271895894526
    6.2982208538220643134764031856321*x1^4-0.84346869683664105199483174146735*x1^3-1.0174703815902903514256649941672*x1^2+0.045382868082479266291784369968809*x1-0.40190925234359098361380802089116*x2^4+2.5405601362437544299410774328862*x2^3-2.5941579088188428947603370033903*x2^2-0.22163699334508613070227056596195*x2+7.9223466274799605457701545674354*x3^4-0.49681209037397611183450862881728*x3^3-4.9218592380105895545128191770345*x3^2-0.39371744715936853042936860447298*x3-0.1971233840446668800217011607856];
gg = [1;1];
input = [0;gg*u1;gg*u2];
sym = [x2-x3^2; x3-x1^2+input(2); -x1-2*x2-x3+x2^3+input(3)];
%%
V = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
% V = 1*x1^4+1*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; C0 = 5.862834287294065;
% V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; 
C0 = 21.873341926982228;
%%
% f = [-x2-3/2*x1^2-1/2*x1^3; x1 - u];
% V = x1^2+x2^2+1*x1*x2;
% C0 = 3.00000001711;
% C0 = 2;
% V = 3*x1^2+2*x2^2+1*x1*x2;
% C0 = 2;
% C0 = 0.1;
%%
% f = [x2-x1
%     0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608
%     ];
% gg = 1;
% V = x1^2+x1*x2+x2^2;
% C0 = 100*8.109410321812767;
%%
dom = 8;
boundary_u = 100;
k = 2;
L_us = 4;
L_au = 4;
gamma = 1;
kk = 1;
i = 0; j = 0;
solh = C0 - V;
trace_Q = 0;
mm = 0;
TRACE = [];
Barrier = [];
%%
C1 = (x1-2)^2+(x2-1)^2+(x3-2)^2-1;
C2 = (x1+1)^2+(x2+2)^2+(x3+1)^2-1;
C3 = (x1+0)^2+(x2-0)^2+(x3-6)^2-4;
C4 = (x1+0)^2+(x2-3)^2+(x3+6)^2-4;
C5 = -x2+2;
C6 = (x1+0.5)^2+(x2+3.25)^2-1;
C7 = -x2 - 4;
C8 = -x1 + 1;
C = [C1;C2;C3;C4;C5;C6;C7;C8];
%%
figure(12);clf;hold on;view(-150, 30);
domain = [-dom dom -dom dom -dom dom];
us1 = patch(pcontour3(C(1),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us1, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us2 = patch(pcontour3(C(2),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us2, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us3 = patch(pcontour3(C(3),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us3, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us4 = patch(pcontour3(C(4),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us4, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
inV = patch(pcontour3(V,C0,domain,'g'));              % Plot the original Lyapunov sublevel set
set(inV, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'g','LineStyle','-','LineWidth',0.7 ); hold on;
%%%%%%%%%%%%%%%%%%%%%%%%
% while abs(double(trace_Q)-double(trace_Q1))>=epsi
% while 1
%
%     mm = mm+1; kk = 1;
%     fprintf('The whole Iteration time is:  %d\n  ',mm);

for i = 1:30
    fprintf('i=%6.0f\n',i);
    if kk == 0
        break
    else
        [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1_3D(f,k,solh,V,mm,gamma,gg,L_au,boundary_u);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if kk == 0
        break
    else
        [solh,trace_Q,kk]=sos_function_2_3D(f,k,SOLu1,SOLu2,SOL1,SOL2,gamma,mm,V,C,dom,gg,L_us);
        TRACE = [TRACE; double(trace_Q)];
        Barrier = [Barrier; solh];
    end
    %%%%%%%
end
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-150, 30);

%     solu =sos_function_3(k,solh,gamma,mm,V,C,V0);
%%
%     axis(domain)
%     kk = 1;
%     while kk == 1
%         j = j + 1;
%         fprintf('j=%6.0f\n',j);
%
%         record = solh;
%         record_Q = trace_Q;
%         [SOLu,SOL1,SOL2] = sos_function_4(k,solh,V,mm,gamma);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [solh, trace_Q, Q, kk]=sos_function_5(k,SOLu,SOL1,SOL2,gamma,mm,V,C);
%         if kk == 0
%             solh = record;
%             trace_Q = record_Q;
%         end
%     %%%%%%%
%     end

%     fprintf('The second round is end.====== \');
% %%
% end