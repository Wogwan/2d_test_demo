clear;
tic
pvar x1 x2 x3 u1 u2;
x = [x1;x2;x3];
%%
% f = [x2-x3^2; x3-x1^2; -x1-2*x2-x3+x2^3];
f = [x3^2+x2
    0.20168850953378745366473268063601*x1^4+0.063529983054143205658580619337485*x1^3-1.5117281788498260486353075293664*x1^2+0.11050308587243914680550528927443*x1-0.010789899544848068069224922282956*x2^4-0.0023351458555991655777206439381644*x2^3-0.050856189012352337464051288407063*x2^2-0.011697431254263418801131457769316*x2-0.0050275698408848868564691159122049*x3^4-0.072414399333034673578168849417125*x3^3+0.043478883248082855761396103844163*x3^2+1.1021325793385576968796968344577*x3-0.41904249821184156449271895894526
    6.2982208538220643134764031856321*x1^4-0.84346869683664105199483174146735*x1^3-1.0174703815902903514256649941672*x1^2+0.045382868082479266291784369968809*x1-0.40190925234359098361380802089116*x2^4+2.5405601362437544299410774328862*x2^3-2.5941579088188428947603370033903*x2^2-0.22163699334508613070227056596195*x2+7.9223466274799605457701545674354*x3^4-0.49681209037397611183450862881728*x3^3-4.9218592380105895545128191770345*x3^2-0.39371744715936853042936860447298*x3-0.1971233840446668800217011607856];
gg = [1;1];
input = [0;gg*u1;gg*u2];
%%
% V = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
V = 6*x1^4-4*x2^4-4*x3^4+4*x1^2*x2^2+2*x1^2*x3^2+4*x2^2*x3^2;
%%
C0 = 1;
cc = 1.2;
k = 4;
dom = 10;
kk = 1;
solU = [];
v_c = [];
iter = 1;
boundary_u = 100;
domain = [-dom dom -dom dom -dom dom];
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
figure(11);clf;hold on;view(-150, 30);
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
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-150, 30);hold on;
%%
while abs(double(cc)-double(C0)) >= 1e-8
    iter = iter + 1;
    if iter ~= 1
        C0 = cc;
    end
    [solu1,solu2,solL,kk]= sos_function_v_3D(f,gg,k,V,C0,boundary_u);
    [cc,kk,solu] = sos_function_v2_3D(f,gg,k,V,C,dom,solL,boundary_u);
    v_c = [v_c; double(cc)];
    solU = [solU;solu];
    figure(11);hold on;
    if kk == 0 && iter == 2
        inV = patch(pcontour3(V,C0,domain,'b'));              % Plot the original Lyapunov sublevel set
        set(inV, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',0.7 ); hold on;
        break
    elseif kk == 0
        inV = patch(pcontour3(V,v_c(end),domain,'b'));              % Plot the original Lyapunov sublevel set
        set(inV, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',0.7 ); hold on;
        break
    end
end
toc

%%
title('');
legend([us4,inV],{'Unsafe Regions','$V_0(x)$'}, 'Interpreter','latex','location','northwest');
xlabel('$x_1$','Interpreter','latex','Fontsize',18);
ylabel('$x_2$','Interpreter','latex','Fontsize',18);
zlabel('$x_3$','Interpreter','latex','Fontsize',18);