%%
figure(10000);clf;hold on;
pvar x1 x2 x3;
dom = 6;
domain = [-dom dom -dom dom -dom dom];
%%
C1 = (x1-3)^2+(x2-0)^2+(x3-0)^2-2;
C2 = (x1+3)^2+(x2+3)^2+(x3-3)^2-3;
C3 = (x1-0)^2+(x2-3)^2+(x3+0)^2-3;
C4 = (x1-3)^2+(x2-0)^2+(x3+3)^2-3;
C = [C1;C2;C3;C4];
us1 = patch(pcontour3(C(1),0,domain,'r'));
set(us1, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us2 = patch(pcontour3(C(2),0,domain,'r'));
set(us2, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us3 = patch(pcontour3(C(3),0,domain,'r'));
set(us3, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us4 = patch(pcontour3(C(4),0,domain,'r'));
set(us4, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
%%
C0 = 0.1;
CC0 = 2.584683714740699;
V = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
inV = patch(pcontour3(V,C0,domain,'b'));
set(inV,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','-','LineWidth',0.7 ); hold on;
%%
inV1 = patch(pcontour3(V,CC0,domain,'b'));
set(inV1,'EdgeAlpha',0.8,'FaceColor', 'none', 'EdgeColor', 'g','LineStyle',':','LineWidth',0.7 ); hold on;
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-150, 30);hold on;
%% New Lyapunov function
V1 = 1.10781433271431*x1^4+0.01711562498101607*x1^3*x2-0.004293891357695184*x1^3*x3-1.819455621599063*x1^2*x2^2-0.07835130081009554*x1^2*x2*x3+0.2005455165134622*x1^2*x3^2+0.005854541235984267*x1*x2^3-0.1053988557975258*x1*x2^2*x3+0.104687646024595*x1*x2*x3^2-0.01576905880117915*x1*x3^3+1.124505629451759*x2^4-0.04474629485339712*x2^3*x3+0.2312786612293422*x2^2*x3^2-0.01730554792698724*x2*x3^3+0.121312765683126*x3^4+0.01156354899908827*x1^3-0.07535038403941373*x1^2*x2-0.01488159063399457*x1^2*x3-0.08501506578247583*x1*x2^2+0.08355970214331113*x1*x2*x3-0.01717793720821389*x1*x3^2+0.01319142674040484*x2^3-0.03173064799350401*x2^2*x3+0.0003880888277869219*x2*x3^2-0.02016579874197291*x3^3+0.3883270529776856*x1^2+0.01952597427976192*x1*x2-0.05128081802475131*x1*x3+0.3766506193975706*x2^2-0.0969329568996595*x2*x3+0.4210891700128911*x3^2+8.835240134064381e-05*x1+0.0003560022249767095*x2-0.0002817478469884745*x3+6.770626500902747;
C1 = 10.296886492629705;
us4 = patch(pcontour3(V1,C1,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us4, 'EdgeAlpha',0.6,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','--','LineWidth',1 ); hold on;

%%
legend([inV,inV1,us4],{'$V(x)\leq c_0$','$V(x)\leq c_1^*$','$V^*(x)\leq c_2^*$'}, 'Interpreter','latex','location','northwest');
title('');
xlabel('$x_1$','Interpreter','latex','Fontsize',18);
ylabel('$x_2$','Interpreter','latex','Fontsize',18);
zlabel('$x_3$','Interpreter','latex','Fontsize',18);
set(gca,'xtick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'ytick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'ztick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'Box','on');view(200,10);
set(gca,'FontSize',22,'Fontname','Times','LineWidth',2);
set(gca,'LooseInset',get(gca,'TightInset'))