%%
figure(16);clf;hold on;
pvar x1 x2 x3;
dom = 10;
domain = [-dom dom -dom dom -dom dom];
%%
C1 = (x1+4)^2+(x2-6)^2+(x3+2)^2-4;
C2 = (x1+3)^2+(x2+4)^2+(x3+4)^2-4;
C3 = (x1-4)^2+(x2-0)^2+(x3-0)^2-5;
C4 = (x1+4)^2+(x2-2)^2+(x3-4)^2-5;
C = [C1;C2;C3;C4];
us1 = patch(pcontour3(C(1),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us1, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us2 = patch(pcontour3(C(2),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us2, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us3 = patch(pcontour3(C(3),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us3, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us4 = patch(pcontour3(C(4),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us4, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;

%%
C0 = 96.811595465770495;
V = 10*x1^4+1*x2^4+20*x3^4+2*x1^2*x2^2-4*x3^2*x2^2+3*x1^2*x3^2;
inV = patch(pcontour3(V,C0,domain,'b'));              % Plot the original Lyapunov sublevel set
set(inV,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','-','LineWidth',1 ); hold on;
%%
B = -11154.51967552443*x1^4+1.135913920002194*x1^3*x2+4.955118594070981*x1^3*x3+1171.840078989759*x1^2*x2^2-15.02990601369866*x1^2*x2*x3+21138.95310790373*x1^2*x3^2-6.189653081344577*x1*x2^3-6.041347173501446*x1*x2^2*x3+6.069007291799152*x1*x2*x3^2+2.099288457169987*x1*x3^3-93.09232391138393*x2^4+3.980296472031688*x2^3*x3-987.2454204626316*x2^2*x3^2+1.819554574488895*x2*x3^3-10085.20085196664*x3^4-3.90264724194604*x1^3-2.546733372542181*x1^2*x2-1.748490559079197*x1^2*x3+1.973311313870531*x1*x2^2+7.025110993477067*x1*x2*x3+8.942932559147787*x1*x3^2+2.947834262954737*x2^3+2.575829927090547*x2^2*x3-1.387124471681909*x2*x3^2-1.835427814430669*x3^3-0.08971485448317769*x1^2+0.3095104865629468*x1*x2+0.2408635828835946*x1*x3-1.1949744015814*x2^2-1.317929814293127*x2*x3-0.4552517850754649*x3^2-0.0001426841292988713*x1+0.0001020831374556339*x2+0.0001228298360055396*x3+998.8215452707433;
B1 = patch(pcontour3(B,0,domain,'b'));              % Plot the original Lyapunov sublevel set
set(B1,'EdgeAlpha',0.6,'FaceColor', 'none', 'EdgeColor', 'g','LineStyle',':','LineWidth',1.6 ); hold on;
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-200, 10);hold on;

%%
legend([inV,B1],{'$B(x)\geq 0$','$B^*(x)\geq 0$'}, 'Interpreter','latex','location','northwest');
title('');
xlabel('$x_1$','Interpreter','latex','Fontsize',18);
ylabel('$x_2$','Interpreter','latex','Fontsize',18);
zlabel('$x_3$','Interpreter','latex','Fontsize',18);
set(gca,'xtick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'ytick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'ztick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'FontSize',22,'Fontname','Times');
set(gca,'Box','on');view(200,10);axis equal;
set(gca,'LooseInset',get(gca,'TightInset'))