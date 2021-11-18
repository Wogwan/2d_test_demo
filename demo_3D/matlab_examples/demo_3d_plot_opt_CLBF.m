%%
figure(15);clf;hold on;
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
V = 0.7754493780412083*x1^4+0.01538143725956307*x1^3*x2-0.01193322783720551*x1^3*x3+0.1406438731320798*x1^2*x2^2+0.1014046006280144*x1^2*x2*x3-0.9550442999604613*x1^2*x3^2+0.0240910069628206*x1*x2^3+0.004303355017772573*x1*x2^2*x3-0.1074291288708021*x1*x2*x3^2-0.04030442091571095*x1*x3^3+0.1044002123665394*x2^4+0.02342017217092511*x2^3*x3+0.3123479952953659*x2^2*x3^2+0.07980855015340288*x2*x3^3+0.6527375651414437*x3^4+0.02009548179139661*x1^3-0.01537810103459676*x1^2*x2-0.01662614893620694*x1^2*x3+0.01738550368743247*x1*x2^2+0.01074722111660453*x1*x2*x3-0.08818121834264446*x1*x3^2+0.001698655876194973*x2^3-0.000237077928820955*x2^2*x3+0.003223829411089607*x2*x3^2+0.009841041529572332*x3^3+0.2675033622283891*x1^2-0.01636121664970925*x1*x2-0.05656214583866465*x1*x3+0.07861645436743456*x2^2+0.0765934355839933*x2*x3+0.1717958491629204*x3^2+0.0004207558255661851*x1-1.769975122612027e-06*x2-4.65115204606246e-05*x3+6.52177325710281;
C0 = 14.961569805334348;
inV = patch(pcontour3(V,C0,domain,'b'));              % Plot the original Lyapunov sublevel set
set(inV,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','-','LineWidth',1 ); hold on;
%%
b1 = -155.7643832928485*x1^4+6.453896364173854*x1^3*x2+3.250614956523293*x1^3*x3+65.71670592337745*x1^2*x2^2+6.188930128566883*x1^2*x2*x3+224.4290112681742*x1^2*x3^2-6.437586351912966*x1*x2^3-5.960365502934631*x1*x2^2*x3+1.282540364805463*x1*x2*x3^2+3.652078784978907*x1*x3^3-21.20939560897872*x2^4+1.01954031910319*x2^3*x3-23.11998086504914*x2^2*x3^2-15.91270159101469*x2*x3^3-100.3009658219718*x3^4-0.4703357597255011*x1^3-2.036667144021659*x1^2*x2-2.479598960735674*x1^2*x3+1.500075456004734*x1*x2^2+7.017102004254787*x1*x2*x3+5.749723534699926*x1*x3^2+2.417892907161381*x2^3+2.311935721992878*x2^2*x3-1.171657570065432*x2*x3^2-0.7141289214602496*x3^3-0.1710268661082545*x1^2+0.1405658844146584*x1*x2+0.153633022415517*x1*x3-0.8165482041800699*x2^2-1.374096821688936*x2*x3-0.7267763655818622*x3^2-0.000441812026141078*x1+0.000287027034801212*x2+0.000340506257482116*x3+1012.112277417705;
B1 = patch(pcontour3(b1,0,domain,'b'));              % Plot the original Lyapunov sublevel set
set(B1,'EdgeAlpha',0.6,'FaceColor', 'none', 'EdgeColor', 'g','LineStyle',':','LineWidth',1.6 ); hold on;
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-200, 10);hold on;

%%
legend([inV,B1],{'$B(x)\geq 0$','$B^*(x)\geq 0$'}, 'Interpreter','latex','location','northwest');
title('');
xlabel('$x_1$','Interpreter','latex','Fontsize',18);
ylabel('$x_2$','Interpreter','latex','Fontsize',18);
zlabel('$x_3$','Interpreter','latex','Fontsize',18);
set(gca,'xtick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'ytick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'ztick',[-dom,-dom/2,0,dom/2,dom]);
set(gca,'FontSize',22,'Fontname','Times');
set(gca,'Box','on');view(200,10);axis equal;
set(gca,'LooseInset',get(gca,'TightInset'))