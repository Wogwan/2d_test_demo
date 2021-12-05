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
sub_level_with_controller = 57.23191152365694;
sublevel_C0 = 2;
%% Fin optimal Lyapunov
optimal_Lya = 0.27761530973907616592910585495702*x1 + 0.28916923872954758412134879108635*x2 - 0.097197146692290778413614305009105*x1^2*x2^2 + 0.16409346533219781871792974925484*x1*x2 - 0.0067107957184185072774251779037513*x1*x2^2 + 0.0039474994197178196742026301535589*x1^2*x2 - 0.014031494938999576269078595203155*x1*x2^3 - 0.00061930885322093666596476868591026*x1^3*x2 + 0.23174538663646407354868017591798*x1^2 - 0.0040444947534192298224664519068483*x1^3 + 0.24289451203744766294434498377086*x2^2 + 0.056072867694129807647485108645924*x1^4 - 0.0014805460738564312894033347944855*x2^3 + 0.066405746174927976488433500890096*x2^4 + 9.4479264108289715551336485077627;
optimal_Lya_sub_level_with_controller = 16.18343934060665;
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
