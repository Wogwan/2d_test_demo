clear;
tic
format long
pvar x1 x2 x3 u1 u2;
x = [x1;x2;x3];
%%
% f = [x2-x3^2; x3-x1^2; -x1-2*x2-x3+x2^3];
f = [x3^2+x2;
    0.030612084516129553501961879646842*x1^4-0.21847098287082121533027363798114*x1^3-0.94331491879401333304896647414266*x1^2+1.0291917885283056825299989138027*x1+0.00055536580209726473921633127517339*x2^4-0.0030111105513248095921774449834629*x2^3-0.002156990191249033225751041698004*x2^2+0.11561511545396291333887006658188*x2-0.13306769730670206519640430542495*x3^4+0.36419745401433389897505321641802*x3^3-0.14211373683420489011375309473806*x3^2+0.9996139813904298291680552979166*x3-0.0046445737253855240262883287586232;
    -0.0045055823774114555496650424970539*x1^4+0.097708541675668683645916701152601*x1^3-x1^2*x3-0.62175706730821389545127431119909*x1^2+0.3968841344326459186220290575875*x1-0.023076088240605339280131502732729*x2^4-0.050675545712074644699729475405547*x2^3+0.047257902328268819314160964495386*x2^2+0.26811627394646608824047007146874*x2-3.780636045475054274334070214536*x3^4+4.4974561702658171213897730922326*x3^3-0.7782608993662247304680026552732*x3^2+0.94957393957775426684975172975101*x3+0.85861129189534040318612003042034];
gg = [1;1];
input = [0;gg*u1;gg*u2];
%%
% V = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
% V = 6*x1^4-4*x2^4-4*x3^4+4*x1^2*x2^2+2*x1^2*x3^2+4*x2^2*x3^2;
V = 4*x1^4-2*x2^4-2*x3^4+1*x1^2*x2^2+6*x1^2*x3^2+5*x2^2*x3^2;
%%
C0 = 1;
cc = 1.2;
k_u = 1;
k_L = 4;
dom = 10;
kk = 1;
solU = [];
v_c = [];
iter = 1;
boundary_u = 200;
domain = [-dom dom -dom dom -dom dom];
%%
C1 = (x1-2)^2+(x2-1)^2+(x3-2)^2-1;
C2 = (x1+1)^2+(x2+2)^2+(x3+1)^2-1;
C3 = (x1+0)^2+(x2-0)^2+(x3-6)^2-4;
C4 = (x1+0)^2+(x2-0)^2+(x3+6)^2-4;
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
    iter = iter + 1
    if iter ~= 1
        C0 = cc;
    end
    [solu1,solu2,solL,kk]= sos_function_v_3D(f,gg,k_u,k_L,V,C0,boundary_u);
    [cc,kk,solu] = sos_function_v2_3D(f,gg,k_u,k_L,V,C,dom,solL,boundary_u);
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