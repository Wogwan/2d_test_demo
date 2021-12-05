%%
figure(10000);clf;hold on;
pvar x1 x2 x3;
dom = 10;
domain = [-dom dom -dom dom -dom dom];
%%
C1 = (x1-3)^2+(x2-0)^2+(x3-0)^2-2;
C2 = (x1+4)^2+(x2+4)^2+(x3-4)^2-4;
C3 = (x1-0)^2+(x2-4)^2+(x3+0)^2-4;
C4 = (x1-4)^2+(x2-0)^2+(x3+4)^2-6;
C = [C1;C2;C3;C4];
% us1 = patch(pcontour3(C(1),0,domain,'r'));
% set(us1, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us2 = patch(pcontour3(C(2),0,domain,'r'));
set(us2, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us3 = patch(pcontour3(C(3),0,domain,'r'));
set(us3, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us4 = patch(pcontour3(C(4),0,domain,'r'));
set(us4, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
%%
C0 = 0.1;
V = x1^4+x1^3*x2+x1^2*x2^2+x1*x2^3+x2^4+x1^3*x3+x1^2*x2*x3+x1*x2^2*x3+x2^3*x3+x1^2*x3^2+x1*x2*x3^2+x2^2*x3^2+x1*x3^3+x2*x3^3+x3^4+x1^3+x1^2*x2+x1*x2^2+x2^3+x1^2*x3+x1*x2*x3+x2^2*x3+x1*x3^2+x2*x3^2+x3^3+x3^2+x1^2+x1*x2+x2^2+x1*x3+x2*x3+x1+x2+x3+1; 
CC0 = 24.92234264283737;inV = patch(pcontour3(V,C0,domain,'b'));
set(inV,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','-','LineWidth',0.7 ); hold on;
%%
inV1 = patch(pcontour3(V,CC0,domain,'b'));
set(inV1,'EdgeAlpha',0.8,'FaceColor', 'none', 'EdgeColor', 'g','LineStyle',':','LineWidth',0.7 ); hold on;
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-150, 30);hold on;
%% New Lyapunov function
V1 = 0.15262216021712840530177857090166*x1 + 0.098262351849115456281502645197179*x2 + 0.12547086987215338993451041460503*x3 + 0.33475379324085863252946637658169*x1^2*x2^2 - 0.33201166533136744485332769727393*x1^2*x3^2 - 1.0336484598631696663062484731199*x2^2*x3^2 + 0.072895564743315710565241261065239*x1*x2 + 0.079779852202941523020562897272612*x1*x3 - 0.028505240380140510481066229431235*x2*x3 + 0.10601990130512847776422802326124*x1*x2^2 - 0.0018975139402669122679762070404763*x1^2*x2 + 0.081849877163025261395112863738177*x1*x2^3 - 0.23231466771492709222357575526985*x1*x3^2 - 0.02811487841558210631909275889484*x1^2*x3 + 0.057337316563754195386515277732542*x1^3*x2 - 0.022553475511912949441617470824895*x1*x3^3 + 0.0024212740083352996336985007985731*x2*x3^2 + 0.029835719308269368493791162677553*x1^3*x3 - 0.075585016094567705757789610743203*x2^2*x3 + 0.057935730725925567441425556580725*x2*x3^3 - 0.047462291878356904772928714919544*x2^3*x3 + 0.095623523343395161466773402025865*x1^2 + 0.038340716855404646801197543481976*x1^3 + 0.10779239865203520121195168712802*x2^2 + 0.11242606670229465803956259151164*x1^4 + 0.022620214673315798548092203645865*x2^3 + 0.049163527809058013562371058924327*x3^2 + 0.45955749986413146901398363297631*x2^4 - 0.0047924211852991144278135493550508*x3^3 + 0.74908482625988148662088406126713*x3^4 - 0.088614880647068336450011827309936*x1*x2*x3^2 + 0.056022378314697314494186031197387*x1*x2^2*x3 - 0.010645205789755573036470259751241*x1^2*x2*x3 + 0.021345023534847486906818048169043*x1*x2*x3 + 20.808884521275363255199408740737;
C1 = 26.55026620001978;
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
set(gca,'Box','on');view(75,8);
set(gca,'FontSize',22,'Fontname','Times','LineWidth',2);
set(gca,'LooseInset',get(gca,'TightInset'))