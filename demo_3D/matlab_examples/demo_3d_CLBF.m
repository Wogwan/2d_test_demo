%%
figure(13);clf;hold on;
pvar x1 x2 x3;
dom = 6;
domain = [-dom dom -dom dom -dom dom];
%%
C1 = (x1+4)^2+(x2-6)^2+(x3+2)^2-4;
C2 = (x1+3)^2+(x2+4)^2+(x3+4)^2-4;
C3 = (x1-4)^2+(x2-0)^2+(x3-0)^2-5;
C4 = (x1+4)^2+(x2-2)^2+(x3-4)^2-5;
C = [C1;C2;C3;C4];
us1 = patch(pcontour3(C(1),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us1, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us2 = patch(pcontour3(C(2),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us2, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us3 = patch(pcontour3(C(3),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us3, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us4 = patch(pcontour3(C(4),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us4, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;

%%
C0 = 13.309586736607509;
V = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
inV = patch(pcontour3(V,C0,domain,'b'));              % Plot the original Lyapunov sublevel set
set(inV,'EdgeAlpha',0.5,'FaceColor', 'none', 'EdgeColor', 'g','LineStyle','-','LineWidth',0.7 ); hold on;
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-150, 30);hold on;

%%
% B = -145.7686630997646*x1^2-18.69634054363406*x1*x2-3.486869598534508*x1*x3-65.61212908891149*x2^2+21.69619390012745*x2*x3-245.3396890062136*x3^2-30.91171455333036*x1+9.877234110946764*x2+121.1541933795845*x3+654.274289233073;
B = (4375004233684943*x2)/9007199254740992 - (5907552702199365*x1)/9007199254740992 + (4901211905026135*x3)/72057594037927936 - (4873286017626133*x1*x2)/4503599627370496 - (1696923333237319*x1*x3)/4503599627370496 - (7321534277155875*x2*x3)/18014398509481984 - (2038086424939521*x1^2)/562949953421312 - (3632306073559465*x2^2)/2251799813685248 - (5362223778392797*x3^2)/1125899906842624 + 5306592666458987/281474976710656;
us4 = patch(pcontour3(B,0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us4, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','--','LineWidth',0.7 ); hold on;

%%
legend([us4,inV],{'$h(x)\geq 0$','$V(x)\leq c_0$'}, 'Interpreter','latex','location','northwest');
title('');
xlabel('$x_1$','Interpreter','latex','Fontsize',18);
ylabel('$x_2$','Interpreter','latex','Fontsize',18);
zlabel('$x_3$','Interpreter','latex','Fontsize',18);
set(gca,'xtick',[-dom,0,dom]);
set(gca,'ytick',[-dom,0,dom]);
set(gca,'ztick',[-dom,0,dom]);
set(gca,'FontSize',18,'Fontname','Times');
set(gca,'Box','on');view(60,10);