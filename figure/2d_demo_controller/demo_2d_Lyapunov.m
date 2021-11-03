% hold on
figure(2);clf;
pvar x1 x2 u htol epsi;
x = [x1;x2];
%%%%%%%%%%%%%%%%%%%%%%%
C1 = (x1+8)^2+(x2-0)^2-4;
C2 = (x1-8)^2+(x2+0)^2-4;
C3 = (x1-0)^2+(x2-8)^2-4;
C4 = (x1-0)^2+(x2+8)^2-4;
V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
sub_level_no_controller = 0.101339366217398;
B_no_controller = -299469.2888551163*x1^6+177574.6581854458*x1^5*x2-779712.3797276674*x1^4*x2^2-868120.9856405794*x1^3*x2^3-1075057.741145297*x1^2*x2^4+524840.3553407069*x1*x2^5-687371.1975566965*x2^6-424386.1054192791*x1^5+2308009.407168391*x1^4*x2-1528276.899024193*x1^3*x2^2-372648.7371172399*x1^2*x2^3-1827051.861052636*x1*x2^4-99505.7212121567*x2^5-1016241.746731024*x1^4+4394543.527266666*x1^3*x2+2156575.298549416*x1^2*x2^2+1064086.444029456*x1*x2^3+180462.8532985262*x2^4-1446177.344812278*x1^3-515132.2066721657*x1^2*x2+4512922.136537937*x1*x2^2+1644329.601763489*x2^3-316821.4412679775*x1^2-5570076.203360715*x1*x2+1441378.478613752*x2^2+962752.0173226857*x1-3245844.443155527*x2+566653.4844229659;
V = 15.22034892819922*x1^4-13.44164154544166*x1^2*x2^2+18.56603884407099*x2^4+5.087438107880891*x1^2+0.2507678441391148*x1*x2+6.951304128544671*x2^2;
H = -12.9254727672629*x1^4+0.02951083237961412*x1^3*x2+10.73667355440355*x1^2*x2^2-1.692619448386333*x1*x2^3-3.582476284173307*x2^4+0.1542660855197659*x1^3+0.06655462590463264*x1^2*x2-0.0764762142639748*x1*x2^2+0.04089981707855826*x2^3+23.41429551310097*x1^2+3.002811313301487*x1*x2-22.87164191870751*x2^2-0.1173283511687331*x1-0.1731079941771329*x2+98.97467757952973;
%%
dom = 14;
domain = [-dom dom -dom dom];
[c24,h24] = pcontour(C1,0,domain,'k'); hold on;    
[~,h21] = pcontour(C2,0,domain,'k'); hold on;    
[~,h22] = pcontour(C3,0,domain,'k'); hold on;            
[~,h23] = pcontour(C4,0,domain,'k'); hold on;       
%%
[c26,h26]=pcontour(V,1.0e+04*1.990871978283739,domain,'r'); hold on; 
[c28,h28]=pcontour(H,0,domain,'g'); hold on; 
[c30,h30]=pcontour(V0,sub_level_no_controller,domain,'m'); hold on; 
[c32,h32]=pcontour(B_no_controller,0,domain,'m'); hold on; 
% clabel(c26,h26);
h24.LineStyle = '--';
h21.LineStyle = '--';
h22.LineStyle = '--';
h23.LineStyle = '--';
h26.LineWidth = 1.2;
h26.LineStyle = '-.';
h26.LineWidth = 1.5;
h28.LineStyle = ':';
h28.LineWidth = 1.8;
h30.LineStyle = '-';
h30.LineWidth = 2.1;
h32.LineStyle = ':';
h32.LineWidth = 1.8;

h = legend([h24,h30,h28,h32,h26],{'Unsafe Regions','$\mathcal{L}_{V_0} = \{0<V(x)\leq 0.1013\}$','$\mathcal{L}_B = \{B(x)\geq 0\}$','$\mathcal{L}_{B_n} = \{B_n(x)\geq 0\}$','$\mathcal{L}_V = \{0<V(x)\leq 19.96\}$'}, 'Interpreter','latex','location','northeast','Fontsize',12);
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal');hold on;
% xlabel('$x_1$','Interpreter','latex','Fontsize',18); 
% ylabel('$x_2$','Interpreter','latex','Fontsize',18);
set(gca,'ytick',[-12,-6,0,6,12]);
set(gca,'xtick',[-12,-6,0,6,12]);
set(gca,'FontSize',14,'Fontname','Times');
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.7;
xlabel('$x_1$','Interpreter','latex','Fontsize',16,'Fontname','Times');
ylabel('$x_2$','Interpreter','latex','Fontsize',16,'Fontname','Times');
xlim([-dom-3 dom+3]); ylim([-dom+2 dom+6]); hold on;

%  Green [54 185 132]/255
%  Red [247 77 77]/255
%  Blue [35 145 213]/255
%  Cyan [121 199 243]/255
%  Grey [133 134 138]/255
