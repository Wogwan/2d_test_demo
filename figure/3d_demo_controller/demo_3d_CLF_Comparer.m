%%
figure(13);clf;hold on;
pvar x1 x2 x3;
dom = 10;
domain = [-dom dom -dom dom -dom dom];
%%
C1 = (x1+4)^2+(x2-6)^2+(x3+2)^2-4;
C2 = (x1+3)^2+(x2+4)^2+(x3+4)^2-4;
C3 = (x1-4)^2+(x2-0)^2+(x3-0)^2-5;
C4 = (x1+4)^2+(x2-2)^2+(x3-4)^2-5;
C = [C1;C2;C3;C4];
us1 = patch(pcontour3(C(1),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us1, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us2 = patch(pcontour3(C(2),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us2, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us3 = patch(pcontour3(C(3),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us3, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;
us4 = patch(pcontour3(C(4),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us4, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','--','LineWidth',0.7 ); hold on;

%%
C0 = 0.1;
C = 96.811595465770495;
V = 10*x1^4+1*x2^4+20*x3^4+2*x1^2*x2^2-4*x3^2*x2^2+3*x1^2*x3^2;
inV = patch(pcontour3(V,C0,domain,'b'));              % Plot the original Lyapunov sublevel set
set(inV,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','-','LineWidth',0.7 ); hold on;
%%
inV1 = patch(pcontour3(V,C,domain,'b'));              % Plot the original Lyapunov sublevel set
set(inV1,'EdgeAlpha',0.8,'FaceColor', 'none', 'EdgeColor', 'g','LineStyle',':','LineWidth',0.7 ); hold on;
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-150, 30);hold on;
%% New Lyapunov function
% V = 0.6208611094349646*x1^4-0.02716156081164378*x1^3*x2-0.02708885812522822*x1^3*x3+0.2240094021882353*x1^2*x2^2+0.4399640327471952*x1^2*x2*x3-0.4691321899483339*x1^2*x3^2+0.04251686032081572*x1*x2^3+0.009187018405654742*x1*x2^2*x3-0.09220713868055416*x1*x2*x3^2-0.06760769468185925*x1*x3^3+0.1488568286343931*x2^4+0.02325223946328157*x2^3*x3+0.3816787511347403*x2^2*x3^2-0.04257663299873744*x2*x3^3+0.467569625242744*x3^4+0.03201239138642915*x1^3-0.01887157059636539*x1^2*x2-0.02573972940739003*x1^2*x3+0.03162825132940131*x1*x2^2+0.04361577201095797*x1*x2*x3-0.1037075490574651*x1*x3^2-0.003366067501592693*x2^3-0.01608336337051472*x2^2*x3+0.01202496046162058*x2*x3^2+0.0233276368120514*x3^3+0.9707911104056866*x1^2-0.01988459016624386*x1*x2-0.1019143611403437*x1*x3+0.3477003513384749*x2^2+0.1331396702460586*x2*x3+0.5228630711587768*x3^2+0.001316609278053213*x1+3.001695459008127e-05*x2-4.3784027128929e-05*x3+9.260142686821981;
% C0 = 18.453381984491568; 
V = 0.7754493780412083*x1^4+0.01538143725956307*x1^3*x2-0.01193322783720551*x1^3*x3+0.1406438731320798*x1^2*x2^2+0.1014046006280144*x1^2*x2*x3-0.9550442999604613*x1^2*x3^2+0.0240910069628206*x1*x2^3+0.004303355017772573*x1*x2^2*x3-0.1074291288708021*x1*x2*x3^2-0.04030442091571095*x1*x3^3+0.1044002123665394*x2^4+0.02342017217092511*x2^3*x3+0.3123479952953659*x2^2*x3^2+0.07980855015340288*x2*x3^3+0.6527375651414437*x3^4+0.02009548179139661*x1^3-0.01537810103459676*x1^2*x2-0.01662614893620694*x1^2*x3+0.01738550368743247*x1*x2^2+0.01074722111660453*x1*x2*x3-0.08818121834264446*x1*x3^2+0.001698655876194973*x2^3-0.000237077928820955*x2^2*x3+0.003223829411089607*x2*x3^2+0.009841041529572332*x3^3+0.2675033622283891*x1^2-0.01636121664970925*x1*x2-0.05656214583866465*x1*x3+0.07861645436743456*x2^2+0.0765934355839933*x2*x3+0.1717958491629204*x3^2+0.0004207558255661851*x1-1.769975122612027e-06*x2-4.65115204606246e-05*x3+6.52177325710281;
C0 = 14.961569805334348;
us4 = patch(pcontour3(V,C0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us4, 'EdgeAlpha',0.6,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','--','LineWidth',1 ); hold on;

%%
legend([inV,inV1,us4],{'$V(x)\leq c_0$','$V(x)\leq c^*$','$V^*(x)\leq c^*$'}, 'Interpreter','latex','location','northwest');
title('');
xlabel('$x_1$','Interpreter','latex','Fontsize',18);
ylabel('$x_2$','Interpreter','latex','Fontsize',18);
zlabel('$x_3$','Interpreter','latex','Fontsize',18);
set(gca,'xtick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'ytick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'ztick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'FontSize',24,'Fontname','Times');
set(gca,'Box','on');view(310,10);axis equal;
set(gca,'LooseInset',get(gca,'TightInset'))