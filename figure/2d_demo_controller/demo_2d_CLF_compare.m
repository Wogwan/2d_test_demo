figure(999);clf;hold on
pcircle(-4,5,2,'k');
pcircle(0,-5,2,'k');
pcircle(5,0,2,'k');
%%
pvar x1 x2 u htol epsi;
x = [x1;x2];
%%
dom1 = 8; domain1 = [-dom1 dom1 -dom1 dom1];
dom2 = 3; domain2 = [-dom2 dom2 -dom2 dom2];
%%%%%%%%%%%%%%%%%%%%%%%
C1 = (x1+4)^2+(x2-5)^2-4;
C2 = (x1+0)^2+(x2+5)^2-4;
C3 = (x1-5)^2+(x2-0)^2-4;
[~,h20] = pcontour(C1,0,domain1,'k'); hold on; h20.LineStyle = '--';     
[~,h21] = pcontour(C2,0,domain1,'k'); hold on; h21.LineStyle = '--';   
[~,h22] = pcontour(C3,0,domain1,'k'); hold on; h22.LineStyle = '--';  
%% With controller
V1 = 1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1;
sub_level_with_controller = 57.23191153176421;
sublevel_C0 = 2;
%% Fin optimal Lyapunov
optimal_Lya = 2.6956787038204343964764575503068*x1 + 2.4744521858387473756124563806225*x2 + 0.2096224907462799214030724215263*x1^2*x2^2 - 0.32267055507215436360723970210529*x1*x2 - 0.082021596385071490753482237323624*x1*x2^2 - 0.1079247740191071058823979456065*x1^2*x2 + 0.17385029105500382495819167161244*x1*x2^3 + 0.13662631844392117419495491503767*x1^3*x2 + 3.0240026272167677134916630166117*x1^2 - 0.025737631683101218349474237356844*x1^3 + 2.7471657151301043242597188509535*x2^2 + 0.11207391066710929716787603638295*x1^4 - 0.054349866126164401991527341806432*x2^3 + 0.12175938649898074284116944454581*x2^4 + 21.036753666558297481969930231571;
optimal_Lya_sub_level_with_controller = 49.60417320108461;
%%
[~,h30]=pcontour(V1,double(sub_level_with_controller),domain1,'m'); hold on; h30.LineStyle = '-.'; h30.LineWidth = 1.4;
[~,h31]=pcontour(V1,double(sublevel_C0),domain1,'m'); hold on; h31.LineStyle = '--'; h31.LineWidth = 1.4;
[~,h33]=pcontour(optimal_Lya,optimal_Lya_sub_level_with_controller,domain1,'m'); h33.LineStyle = '-'; h33.LineWidth = 1.8;
%%
h = legend([h22,h31,h30,h33],{'Unsafe Regions','CLF $c_0$','CLF $c^*$','Optimal CLF'}, 'Interpreter','latex','location','northeast','Fontsize',12);
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal');hold on;
set(gca,'ytick',[-10,-5,0,5,10]);
set(gca,'xtick',[-10,-5,0,5,10]);
set(gca,'FontSize',14,'Fontname','Times');
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.7;
xlabel('$x_1$','Interpreter','latex','Fontsize',18,'Fontname','Times');
ylabel('$x_2$','Interpreter','latex','Fontsize',18,'Fontname','Times');
xlim([-dom1+1 dom1+1]); ylim([-dom1 dom1+1]); hold on;
set(gca,'LooseInset',get(gca,'TightInset'))
%%
%  Green [54 185 132]/255
%  Red [247 77 77]/255
%  Blue [35 145 213]/255
%  Cyan [121 199 243]/255
%  Grey [133 134 138]/255
