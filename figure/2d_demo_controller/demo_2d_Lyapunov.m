% hold on
figure(2);clf;
pvar x1 x2 u htol epsi;
x = [x1;x2];
%%%%%%%%%%%%%%%%%%%%%%%
C1 = (x1-0)^2+(x2-2.8)^2-1;
C2 = (x1+3)^2+(x2+2)^2-1;
C3 = (x1-3)^2+(x2+0)^2-1;
C4 = (x1-2)^2+(x2-6)^2-1;
V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
H = -12.9254727672629*x1^4+0.02951083237961412*x1^3*x2+10.73667355440355*x1^2*x2^2-1.692619448386333*x1*x2^3-3.582476284173307*x2^4+0.1542660855197659*x1^3+0.06655462590463264*x1^2*x2-0.0764762142639748*x1*x2^2+0.04089981707855826*x2^3+23.41429551310097*x1^2+3.002811313301487*x1*x2-22.87164191870751*x2^2-0.1173283511687331*x1-0.1731079941771329*x2+98.97467757952973;
%%
dom = 5;
domain = [-dom dom -dom dom];
xlim([-dom dom]); ylim([-dom+1 dom+2]); hold on;
[c24,h24] = pcontour(C1,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
[~,~] = pcontour(C2,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
[~,~] = pcontour(C3,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
[c26,h26]=pcontour(V,19.962969754795161,domain,'r'); hold on; 
[c28,h28]=pcontour(H,0,domain,'g'); hold on; 
% clabel(c26,h26);
h24.LineStyle = '-';
h26.LineWidth = 1.2;
h26.LineStyle = '--';
h26.LineWidth = 1.5;
h28.LineStyle = '-.';
h28.LineWidth = 1.8;

h = legend([h24,h26,h28],{'Unsafe Regions','$\mathcal{L}_V = \{0<V(x)\leq 19.96\}$','$\mathcal{L}_B = \{B(x)\geq 0\}$'}, 'Interpreter','latex','location','northeast','Fontsize',26);
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal');hold on;
% xlabel('$x_1$','Interpreter','latex','Fontsize',18); 
% ylabel('$x_2$','Interpreter','latex','Fontsize',18);
set(gca,'ytick',[-5,-3,0,3,5]);
set(gca,'xtick',[-5,-3,0,3,5]);
set(gca,'FontSize',20,'Fontname','Times');
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.7;
xlabel('$x_1$','Interpreter','latex','Fontsize',22,'Fontname','Times');
ylabel('$x_2$','Interpreter','latex','Fontsize',22,'Fontname','Times');

%  Green [54 185 132]/255
%  Red [247 77 77]/255
%  Blue [35 145 213]/255
%  Cyan [121 199 243]/255
%  Grey [133 134 138]/255
