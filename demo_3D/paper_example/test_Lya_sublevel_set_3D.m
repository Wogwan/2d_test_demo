clear;
tic
pvar x1 x2 x3 u1 u2;
x = [x1;x2;x3];
%%
% f = [x2-x3^2; x3-x1^2; -x1-2*x2-x3+x2^3];
f = [x3^2+x2;
    -0.0049547057928618194316382937132619*x1^6+0.016038104182118455834786776862194*x1^5+0.17304421449199596376975923763518*x1^4-0.015148831221695625174699448385663*x1^3-1.9664452555228579348805786588785*x1^2-0.0021576962449399269430721510408662*x1-0.000011706205756031373211353135976864*x2^6-0.00023592386382095546413194264712132*x2^5+0.0013377477474279379256183464264041*x2^4+0.020916317493555061646226533866866*x2^3+0.0035615524584319516350483514344205*x2^2+0.0024604503765848799098914234662061*x2-0.000050670167990751247867609041719561*x3^6-0.000058840446465146536113299813308686*x3^5-0.000055134099666753736313125344725705*x3^4-0.000033601223035472935836594221559182*x3^3+0.000020495395474796491964118716477827*x3^2+1.000199970340097577106014403725*x3-0.0067744531592688922872619317061549;
    -0.00021427187993059508531144830012494*x1^6-0.0035301243133028783680038564796178*x1^5-0.011827602656307062076179725806924*x1^4+0.011372844865827961419180169855281*x1^3-x1^2*x3-0.0065379100094886729438448114137827*x1^2+0.0016946159218526283141842414536882*x1-0.000055678548612553848101269088344267*x2^6+0.00018812390950980931820721298031174*x2^5+0.0032182348530724823841564496973433*x2^4-0.0085648497718256149519033826322811*x2^3-0.0040248184553112534289631696537981*x2^2-0.003972827938322501251100504759961*x2+18.489742019095768758157882771798*x3^6-2.953382702585764876986040849971*x3^5-31.887235740810809992221153130743*x3^4+1.3823115018389194883288717541348*x3^3+12.624803455757429755841603213895*x3^2-1.0204265342901299336297871178658*x3-0.50995581847337662149792503196721];
gg = [1;1];
input = [0;gg*u1;gg*u2];
%%
% V = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
% V = 6*x1^4-4*x2^4-4*x3^4+4*x1^2*x2^2+2*x1^2*x3^2+4*x2^2*x3^2;
V = 1*x1^4+4*x2^4+1*x3^4+2*x1^2*x2^2+2*x1^2*x3^2+2*x2^2*x3^2;
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