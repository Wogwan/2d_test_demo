clear;
tic
format long
pvar x1 x2 x3 u1 u2;
x = [x1;x2;x3];
%%
% f = [x2-x3^2; x3-x1^2; -x1-2*x2-x3+x2^3];
f = [x3^2+x2;
    0.032985569548582671650649444927694*x1^4-0.25190811831004457417622406696696*x1^3-0.87515455184313891217181756710402*x1^2+0.98308426849193212652503840824162*x1+0.0064285707919207223451363297783701*x2^4+0.0078332491243005321346348779343316*x2^3-0.010712089673140650150640063031915*x2^2+0.053882908343556683294917064586116*x2-0.12661718463574700432872077726643*x3^4+0.29574073949289136908902264622157*x3^3-0.0045483349953734054579856938005378*x3^2+0.87242974766040720657755969114078*x3-0.018123100035730797166089643002492;
    -0.0045563969102914974457219088321835*x1^4+0.098282761847104382901818553364137*x1^3-x1^2*x3-0.63059661735708949503731446384336*x1^2+0.39788196631989014573349550119019*x1-0.022706448443695230465788625906498*x2^4-0.05009258712024410725716094816562*x2^3+0.042410939116867982234815315223386*x2^2+0.26685438332235023040084342937917*x2-3.7823648546322456986956694890978*x3^4+4.5017165904882983085144587676041*x3^3-0.78859871736978876704104137550797*x3^2+0.95794780260139500427385428338312*x3+0.91223086905465527163918611344087];
gg = [1;1];
input = [0;gg*u1;gg*u2];
%%
V = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
% V = 6*x1^4-4*x2^4-4*x3^4+4*x1^2*x2^2+2*x1^2*x3^2+4*x2^2*x3^2;
% V = 4*x1^4-2*x2^4-2*x3^4+1*x1^2*x2^2+6*x1^2*x3^2+5*x2^2*x3^2;
% V = 2*x1^4+2*x2^4+1*x3^4+4*x1^2*x2^2+4*x1^2*x3^2+8*x2^2*x3^2;
%%
C0 = 1;
cc = 1.2;
k_u = 2;
k_L = 4;
dom = 6;
kk = 1;
solU = [];
v_c = [];
iter = 1;
boundary_u = 50;
domain = [-dom dom -dom dom -dom dom];
%%
C1 = (x1-3)^2+(x2-3)^2+(x3-0)^2-3;
C2 = (x1+3)^2+(x2+3)^2+(x3-0)^2-3;
C3 = (x1+0)^2+(x2-0)^2+(x3-4)^2-4;
C4 = (x1+0)^2+(x2-0)^2+(x3+4)^2-4;
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
while abs(double(cc)-double(C0)) >= 1e-4
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