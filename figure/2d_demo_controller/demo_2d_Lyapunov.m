% hold on
figure(2);clf;
pvar x1 x2 u htol epsi;
x = [x1;x2];
%%
dom1 = 10; domain1 = [-dom1 dom1 -dom1 dom1];
dom2 = 3; domain2 = [-dom2 dom2 -dom2 dom2];
%%%%%%%%%%%%%%%%%%%%%%%
C1 = (x1+4)^2+(x2-6)^2-4;
C2 = (x1+3)^2+(x2+4)^2-4;
C3 = (x1-6)^2+(x2-0)^2-5;
[~,h20] = pcontour(C1,0,domain1,'k'); hold on; h20.LineStyle = '--';     
[~,h21] = pcontour(C2,0,domain1,'k'); hold on; h21.LineStyle = '--';   
[~,h22] = pcontour(C3,0,domain1,'k'); hold on; h22.LineStyle = '--';  
% pcircle(-4,6,2,[133 134 138]/255);
% pcircle(-3,-4,2,[133 134 138]/255);
% pcircle(6,0,sqrt(5),[133 134 138]/255);
%% Without controller
V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
sub_level_no_controller = 1e-2;
sub_maximum = 0.3669112237857279;
B_no_controller = -0.01730994859216674*x1^4+0.01942272684429362*x1^3*x2+0.03417833273172551*x1^2*x2^2-0.03089576272906369*x1*x2^3-0.03229854756848641*x2^4+0.0469049245134507*x1^3-0.02724130274950444*x1^2*x2-0.1920732549990291*x1*x2^2-0.04355589754110674*x2^3-0.1908212359824433*x1^2+0.04850161666128519*x1*x2-0.3902288519748806*x2^2+0.02342000091755831*x1+0.01847082505198164*x2+0.1049276784359924;
% [~,h26]=pcontour(V0,double(sub_level_no_controller),domain2,'r'); hold on; h26.LineWidth = 1.5; % h26.LineStyle = '-'; 
% [~,h27]=pcontour(V0,double(sub_maximum),domain2,'r'); hold on; h27.LineWidth = 1.5; % h27.LineStyle = '-'; 
[~,h28]=pcontour(B_no_controller,0,domain2,'r'); hold on; h28.LineStyle = '-'; h28.LineWidth = 1.4;

%% With controller
V1 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
initial_set = 0.1;
sub_level_with_controller = 1e+2*1.237395136833918;
B_controller_1 = -25.5286775994005*x1^4-8.826737822184512*x1^3*x2+9.97662491998584*x1^2*x2^2-57.03373052165595*x1*x2^3-66.6347601541668*x2^4+6.719666264862423*x1^3+27.82600314594399*x1^2*x2+63.24796692083677*x1*x2^2+32.27190868643367*x2^3-12.33737345606981*x1^2-16.10987157050511*x1*x2-5.888600774219452*x2^2-0.003646379657348746*x1-0.002562037970909413*x2+4830.062776216545;
%% Fin optimal Lyapunov
optimal_Lya = 0.04179276814237692*x1^4+0.00599499247615957*x1^3*x2-0.01159727730969628*x1^2*x2^2+0.1133535410885831*x1*x2^3+0.1278491306857938*x2^4-0.02172406709266159*x1^3-0.01538819184502436*x1^2*x2-0.0696562899336908*x1*x2^2-0.06418623865759861*x2^3+1.008634051772913*x1^2-0.1426054475471644*x1*x2+1.639478142955251*x2^2+0.0001459131173134332*x1+0.0002096383926257461*x2+16.81346393865286;
optimal_Lya_initial_sub_level_with_controller = 16.913463938652860;
optimal_Lya_sub_level_with_controller = 36.684617121664921;
B_controller_2 = -0.1905010593814771*x1^4+0.03447036469657989*x1^3*x2+0.262456315960618*x1^2*x2^2-0.2777325963053551*x1*x2^3-0.3466812824113252*x2^4+0.04686504039624215*x1^3+0.3998864275022916*x1^2*x2+1.077592296654437*x1*x2^2+0.3939104597524347*x2^3-0.6938361093119657*x1^2-0.7848739465061079*x1*x2-1.191985495239705*x2^2-0.000218085276651845*x1-0.0003126335767281804*x2+44.08881714176047;
%%
% [~,h29]=pcontour(V1,double(initial_set),domain1,'r'); hold on; h29.LineStyle = '-'; h29.LineWidth = 2.1;
[~,h30]=pcontour(V1,double(sub_level_with_controller),domain1,'m'); hold on; h30.LineStyle = '-.'; h30.LineWidth = 1.4;
[~,h31]=pcontour(B_controller_1,0,domain1,'g'); hold on; h31.LineStyle = '-.'; h31.LineWidth = 2.5;
% [~,h32]=pcontour(optimal_Lya,optimal_Lya_initial_sub_level_with_controller,domain1,'b'); h32.LineStyle = '-'; h32.LineWidth = 2.1;
[~,h33]=pcontour(optimal_Lya,optimal_Lya_sub_level_with_controller,domain1,'m'); h33.LineStyle = '-'; h33.LineWidth = 1.8;
[~,h34]=pcontour(B_controller_2,0,domain1,'b'); h34.LineStyle = '-'; h34.LineWidth = 2.5;
%%
h = legend([h22,h28,h30,h33,h31,h34],{'Unsafe Regions','BC','CLF','Optimal CLF','CBF','Optimal CBF'}, 'Interpreter','latex','location','northeast','Fontsize',12);
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal');hold on;
% xlabel('$x_1$','Interpreter','latex','Fontsize',18); 
% ylabel('$x_2$','Interpreter','latex','Fontsize',18);
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

%  Green [54 185 132]/255
%  Red [247 77 77]/255
%  Blue [35 145 213]/255
%  Cyan [121 199 243]/255
%  Grey [133 134 138]/255
