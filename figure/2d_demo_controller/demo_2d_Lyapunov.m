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
sub_level_with_controller = 1.0e+02*1.237395136880424;
B_controller_1 = 9.5200727248137635427838176838122*x1^2*x2^2-0.0019067415715053864094796765016326*x2-0.0027124382181331805306834237256908*x1-16.59753686623925972298820852302*x1*x2+63.590614823783475628715677885339*x1*x2^2+27.736623424839930152074884972535*x1^2*x2-57.175308245138516838323994306847*x1*x2^3-8.9090507766416919821494957432151*x1^3*x2-12.500188005127661483584233792499*x1^2+6.8887349167723126441842396161519*x1^3-6.1066090718013157356836018152535*x2^2-25.192936328667038026196678401902*x1^4+32.950912792629338810002082027495*x2^3-66.857596518457086176567827351391*x2^4+4858.2399374571323278360068798065;
%% Fin optimal Lyapunov
optimal_Lya = 0.04347729711111809*x1^4+0.002521216040839803*x1^3*x2-0.009578250553215793*x1^2*x2^2+0.1107974772000051*x1*x2^3+0.1326098154646251*x2^4-0.02107874625877001*x1^3-0.00520607188316635*x1^2*x2-0.06500720514263179*x1*x2^2-0.05111584728043024*x2^3+1.242392502164554*x1^2-0.03874961976499153*x1*x2+1.999865552141553*x2^2+0.0001861635321948478*x1+0.0002003093521683895*x2+18.28135511117815;
optimal_Lya_sub_level_with_controller = 41.217845170661541;
B_controller_2 = -0.2986730272815503*x1^4+0.03683820773584286*x1^3*x2+0.5163488949428232*x1^2*x2^2-0.2518904549547729*x1*x2^3-0.451106485172958*x2^4+0.1431615752994847*x1^3+0.3396504971750378*x1^2*x2+0.9681904880887063*x1*x2^2+0.493339341250861*x2^3-0.8568200107939279*x1^2-0.9026805442363007*x1*x2-1.682519849229286*x2^2-0.0001618032041516801*x1-0.0002318316938755713*x2+46.45147306706053;
%%
[~,h30]=pcontour(V1,double(sub_level_with_controller),domain1,'m'); hold on; h30.LineStyle = '-.'; h30.LineWidth = 1.4;
[~,h31]=pcontour(B_controller_1,0,domain1,'g'); hold on; h31.LineStyle = '-.'; h31.LineWidth = 2.5;
[~,h33]=pcontour(optimal_Lya,optimal_Lya_sub_level_with_controller,domain1,'m'); h33.LineStyle = '-'; h33.LineWidth = 1.8;
[~,h34]=pcontour(B_controller_2,0,domain1,'b'); h34.LineStyle = '-'; h34.LineWidth = 2.5;
%%
h = legend([h22,h28,h30,h33,h31,h34],{'Unsafe Regions','BC','CLF','Optimal CLF','CBF','Optimal CBF'}, 'Interpreter','latex','location','northeast','Fontsize',12);
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
<<<<<<< HEAD
xlim([-dom1+1 dom1+3]); ylim([-dom1 dom1+2]); hold on;
set(gca,'LooseInset',get(gca,'TightInset'))
=======
xlim([-dom1+2 dom1+3]); ylim([-dom1 dom1+7.5]); hold on;

>>>>>>> 4c7f0ecb8a2f24b865afd7b40c7f24b87f92ad16
%  Green [54 185 132]/255
%  Red [247 77 77]/255
%  Blue [35 145 213]/255
%  Cyan [121 199 243]/255
%  Grey [133 134 138]/255
