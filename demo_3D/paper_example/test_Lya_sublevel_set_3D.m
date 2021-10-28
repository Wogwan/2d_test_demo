clear;
tic
pvar x1 x2 x3 u1 u2;
x = [x1;x2;x3];
%%
% f = [x2-x3^2; x3-x1^2; -x1-2*x2-x3+x2^3];
f = [x3^2+x2;
    -0.013439783370636955439625381814039*x1^6-0.016442956217049792960738230362949*x1^5+0.03935745159850712136112571570834*x1^4+0.023551174673843555154072639306363*x1^3-1.416309426151469999521914644447*x1^2+0.038264882251191960335493974833582*x1+0.0021070037035007486286852795842606*x2^6+0.010008081632901627208709349758919*x2^5+0.021392236846794521892833884635365*x2^4-0.0084261693046387264177665699094177*x2^3-0.010830198130138990811333066233146*x2^2+0.014805390770193117486175360397738*x2+0.0018954677780157298209312566328322*x3^6-0.010500768869064850200012450898157*x3^5-0.0017268337014203910539933417567227*x3^4+0.079777984342271290874037958928966*x3^3-0.045308072920140940453848088509403*x3^2+1.0084310754092052873909235444216*x3-0.47055648698196069976140698543077;
    0.37622655637946228468493359287095*x1^6-0.06211212667655827135426704899146*x1^5-0.052163923784865229293927768594585*x1^4+0.081140591279098853161322324467619*x1^3-0.077950826976361506370771792262531*x1^2-0.93405193717951635890006656381956*x1-0.020276638443425185065471794132463*x2^6-0.026242945812779666647784893029893*x2^5-0.12538335329298685993926198989357*x2^4+1.1492783526438503927113998770437*x2^3-0.11154140442300677915632434178406*x2^2-1.9297474994111820240094701262024*x2+0.81796822550108683191893987896037*x3^6+0.068966467071334303096108442332479*x3^5+1.102306576917776531621129265659*x3^4-4.4400325167154852390449804033778*x3^3-4.0498079986114728645585358890457*x3^2-0.33174966696044531910825270415444*x3-0.37205321834097124927831501395303];
gg = [1;1];
input = [0;gg*u1;gg*u2];
%%
% V = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+2*x1^2*x3^2+2*x2^2*x3^2;
V = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
% V = 3*x1^2+10*x2^2+10*x3^2+12*x1*x2+2*x1*x3+6*x2*x3;
%%
C0 = 1;
cc = 1.2;
k = 6;
dom = 10;
kk = 1;
solU = [];
v_c = [];
iter = 1;
boundary_u = 100;
domain = [-dom dom -dom dom -dom dom];
%%
C1 = (x1-2)^2+(x2-1)^2+(x3-2)^2-1;
C2 = (x1+1)^2+(x2+2)^2+(x3+1)^2-1;
C3 = (x1+0)^2+(x2-0)^2+(x3-6)^2-4;
C4 = (x1+0)^2+(x2-3)^2+(x3+6)^2-4;
C5 = -x2+2;
C6 = (x1+0.5)^2+(x2+3.25)^2-1;
C7 = -x2 - 4;
C8 = -x1 + 1;
C = [C1;C2;C3;C4;C5;C6;C7;C8];
%%
figure(11);clf;hold on;view(-150, 30);
us1 = patch(pcontour3(C(1),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us1, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us2 = patch(pcontour3(C(2),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us2, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us3 = patch(pcontour3(C(3),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us3, 'EdgeAlpha',0.05,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
us4 = patch(pcontour3(C(4),0,domain,'r'));            % Plot the original Lyapunov sublevel set
set(us4, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'r','LineStyle','--','LineWidth',0.7 ); hold on;
inV = patch(pcontour3(V,C0,domain,'g'));              % Plot the original Lyapunov sublevel set
set(inV, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'g','LineStyle','-','LineWidth',0.7 ); hold on;
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-150, 30);hold on; 
%%
while abs(double(cc)-double(C0)) >= 1e-8
    iter = iter + 1;
    if iter ~= 1
        C0 = cc;
    end
    [solu1,solu2,solL,kk]= sos_function_v_3D(f,gg,k,V,C0,boundary_u);
    [cc,kk,solu] = sos_function_v2_3D(f,gg,k,V,C,dom,solL,boundary_u);
    v_c = [v_c; double(cc)];
    solU = [solU;solu];
    figure(11);hold on;
    if kk == 0 && iter == 2
        inV = patch(pcontour3(V,C0,domain,'b'));              % Plot the original Lyapunov sublevel set
        set(inV, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',0.7 ); hold on;
        break
    elseif kk == 0
        inV = patch(pcontour3(V,v_c(end),domain,'b'));              % Plot the original Lyapunov sublevel set
        set(inV, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',0.7 ); hold on;
        break
    end
end
toc

%%
title('');
legend([us4,inV],{'Unsafe Regions','$V_0(x)$'}, 'Interpreter','latex','location','northwest');
xlabel('$x_1$','Interpreter','latex','Fontsize',18);
ylabel('$x_2$','Interpreter','latex','Fontsize',18);
zlabel('$x_3$','Interpreter','latex','Fontsize',18);