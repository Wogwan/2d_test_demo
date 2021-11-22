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

%% Threshold 1e-9 Version
% V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; %% Without controller
% sub_level_no_controller = 1e-2;
% sub_maximum = 0.3669112237857279;
% B_no_controller = -0.01730994859216674*x1^4+0.01942272684429362*x1^3*x2+0.03417833273172551*x1^2*x2^2-0.03089576272906369*x1*x2^3-0.03229854756848641*x2^4+0.0469049245134507*x1^3-0.02724130274950444*x1^2*x2-0.1920732549990291*x1*x2^2-0.04355589754110674*x2^3-0.1908212359824433*x1^2+0.04850161666128519*x1*x2-0.3902288519748806*x2^2+0.02342000091755831*x1+0.01847082505198164*x2+0.1049276784359924;
% % [~,h26]=pcontour(V0,double(sub_level_no_controller),domain2,'r'); hold on; h26.LineWidth = 1.5; % h26.LineStyle = '-'; [~,h27]=pcontour(V0,double(sub_maximum),domain2,'r'); hold on; h27.LineWidth = 1.5; % h27.LineStyle = '-'; 
% [~,h28]=pcontour(B_no_controller,0,domain2,'r'); hold on; h28.LineStyle = '-'; h28.LineWidth = 1.4;
% V1 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; %% With controller
% initial_set = 0.1;
% sub_level_with_controller = 1.0e+02*1.237395136812358;
% B_controller_1 = -25.5286775994005*x1^4-8.826737822184512*x1^3*x2+9.97662491998584*x1^2*x2^2-57.03373052165595*x1*x2^3-66.6347601541668*x2^4+6.719666264862423*x1^3+27.82600314594399*x1^2*x2+63.24796692083677*x1*x2^2+32.27190868643367*x2^3-12.33737345606981*x1^2-16.10987157050511*x1*x2-5.888600774219452*x2^2-0.003646379657348746*x1-0.002562037970909413*x2+4830.062776216545;

%% Threshold 1e-6 Version
V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; %% Without controller
sub_level_no_controller = 1e-2;
sub_maximum = 0.3669112237857279;
B_no_controller = -0.01730994859216674*x1^4+0.01942272684429362*x1^3*x2+0.03417833273172551*x1^2*x2^2-0.03089576272906369*x1*x2^3-0.03229854756848641*x2^4+0.0469049245134507*x1^3-0.02724130274950444*x1^2*x2-0.1920732549990291*x1*x2^2-0.04355589754110674*x2^3-0.1908212359824433*x1^2+0.04850161666128519*x1*x2-0.3902288519748806*x2^2+0.02342000091755831*x1+0.01847082505198164*x2+0.1049276784359924;
% [~,h26]=pcontour(V0,double(sub_level_no_controller),domain2,'r'); hold on; h26.LineWidth = 1.5; % h26.LineStyle = '-'; [~,h27]=pcontour(V0,double(sub_maximum),domain2,'r'); hold on; h27.LineWidth = 1.5; % h27.LineStyle = '-'; 
[~,h28]=pcontour(B_no_controller,0,domain2,'r'); hold on; h28.LineStyle = '-'; h28.LineWidth = 1.4;
%% With controller
V1 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
initial_set = 0.1;
sub_level_with_controller = 1.0e+02*1.237395136880424;
B_controller_1 = -25.19293632866704*x1^4-8.909050776641692*x1^3*x2+9.520072724813764*x1^2*x2^2-57.17530824513852*x1*x2^3-66.85759651845709*x2^4+6.888734916772313*x1^3+27.73662342483993*x1^2*x2+63.59061482378348*x1*x2^2+32.95091279262934*x2^3-12.50018800512766*x1^2-16.59753686623926*x1*x2-6.106609071801316*x2^2-0.002712438218133181*x1-0.001906741571505386*x2+4858.239937457132;
%% Fin optimal Lyapunov
optimal_Lya = 0.0492406555081527*x1^4+0.0002237009592458295*x1^3*x2-0.01435222949115033*x1^2*x2^2+0.1218754852041238*x1*x2^3+0.148629315170593*x2^4-0.02511065097336655*x1^3-0.004003280470573109*x1^2*x2-0.06608123542663244*x1*x2^2-0.05056025999262872*x2^3+1.253847299038444*x1^2+0.06087863420318874*x1*x2+1.88511390003792*x2^2+0.0002352778273053929*x1+0.0002504525971563049*x2+7.604033418627425;
optimal_Lya_sub_level_with_controller = 31.131657156172718; 
B_controller_2 = -0.2528928353460989*x1^4+0.08243074581745588*x1^3*x2+0.4217151522820832*x1^2*x2^2-0.2971362743864769*x1*x2^3-0.4218913155007792*x2^4+0.1617303358687888*x1^3+0.4712288380128266*x1^2*x2+1.143485282202283*x1*x2^2+0.4275270563246917*x2^3-0.8243376278612055*x1^2-1.127048576999814*x1*x2-1.543607394766479*x2^2-0.0002390944157163832*x1-0.0003378347453100939*x2+49.81717543510401;

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
set(gca,'LooseInset',get(gca,'TightInset'))
%  Green [54 185 132]/255
%  Red [247 77 77]/255
%  Blue [35 145 213]/255
%  Cyan [121 199 243]/255
%  Grey [133 134 138]/255
