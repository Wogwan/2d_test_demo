figure(999);clf;hold on
pcircle(-4,6,2,'k');
pcircle(-3,-4,2,'k');
pcircle(6,0,sqrt(5),'k');
%%
pvar x1 x2 u htol epsi;
x = [x1;x2];
%%
dom1 = 8; domain1 = [-dom1 dom1 -dom1 dom1];
dom2 = 3; domain2 = [-dom2 dom2 -dom2 dom2];
%%%%%%%%%%%%%%%%%%%%%%%
C1 = (x1+4)^2+(x2-6)^2-4;
C2 = (x1+3)^2+(x2+4)^2-4;
C3 = (x1-6)^2+(x2-0)^2-5;
[~,h20] = pcontour(C1,0,domain1,'k'); hold on; h20.LineStyle = '--';     
[~,h21] = pcontour(C2,0,domain1,'k'); hold on; h21.LineStyle = '--';   
[~,h22] = pcontour(C3,0,domain1,'k'); hold on; h22.LineStyle = '--';  
%% Without Controller
B_no_controller = -0.01730994859216674*x1^4+0.01942272684429362*x1^3*x2+0.03417833273172551*x1^2*x2^2-0.03089576272906369*x1*x2^3-0.03229854756848641*x2^4+0.0469049245134507*x1^3-0.02724130274950444*x1^2*x2-0.1920732549990291*x1*x2^2-0.04355589754110674*x2^3-0.1908212359824433*x1^2+0.04850161666128519*x1*x2-0.3902288519748806*x2^2+0.02342000091755831*x1+0.01847082505198164*x2+0.1049276784359924;
[~,h28]=pcontour(B_no_controller,0,domain2,'r'); hold on; h28.LineStyle = '-'; h28.LineWidth = 1.4;
%% With controller
V1 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
sub_level_with_controller = 89.952634023568450;
%% Fin optimal Lyapunov
optimal_Lya = 0.04375278098272024*x1^4+0.002929993727326272*x1^3*x2-0.009153369646952632*x1^2*x2^2+0.1117128919832005*x1*x2^3+0.1337597502643151*x2^4-0.02068552493110127*x1^3-0.005144335364106218*x1^2*x2-0.06593999449113327*x1*x2^2-0.05150085197129949*x2^3+1.270472041414144*x1^2-0.04520082154572778*x1*x2+2.036327882785228*x2^2+0.0001597265387810012*x1+0.0001773175441663645*x2+18.32229824450674;
optimal_Lya_sub_level_with_controller = 41.61381252856700;
%%
[~,h30]=pcontour(V1,double(sub_level_with_controller),domain1,'m'); hold on; h30.LineStyle = '-.'; h30.LineWidth = 1.4;
[~,h33]=pcontour(optimal_Lya,optimal_Lya_sub_level_with_controller,domain1,'m'); h33.LineStyle = '-'; h33.LineWidth = 1.8;
%%
h = legend([h22,h30,h33],{'Unsafe Regions','CLF','Optimal CLF'}, 'Interpreter','latex','location','northeast','Fontsize',12);
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
xlim([-dom1+1 dom1+3]); ylim([-dom1 dom1+2]); hold on;
set(gca,'LooseInset',get(gca,'TightInset'))
%%
%  Green [54 185 132]/255
%  Red [247 77 77]/255
%  Blue [35 145 213]/255
%  Cyan [121 199 243]/255
%  Grey [133 134 138]/255
