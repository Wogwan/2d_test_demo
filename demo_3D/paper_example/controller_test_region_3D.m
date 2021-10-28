clear;
pvar x1 x2 x3 u1 u2 htol epsi;
x = [x1;x2;x3];
%%%%%%%%%%%%%%%%%%%%%%%
%%
% f = [x2-x3^2; x3-x1^2; -x1-2*x2-x3+x2^3];
f = [x3^2+x2;
    -0.013439783370636955439625381814039*x1^6-0.016442956217049792960738230362949*x1^5+0.03935745159850712136112571570834*x1^4+0.023551174673843555154072639306363*x1^3-1.416309426151469999521914644447*x1^2+0.038264882251191960335493974833582*x1+0.0021070037035007486286852795842606*x2^6+0.010008081632901627208709349758919*x2^5+0.021392236846794521892833884635365*x2^4-0.0084261693046387264177665699094177*x2^3-0.010830198130138990811333066233146*x2^2+0.014805390770193117486175360397738*x2+0.0018954677780157298209312566328322*x3^6-0.010500768869064850200012450898157*x3^5-0.0017268337014203910539933417567227*x3^4+0.079777984342271290874037958928966*x3^3-0.045308072920140940453848088509403*x3^2+1.0084310754092052873909235444216*x3-0.47055648698196069976140698543077;
    0.37622655637946228468493359287095*x1^6-0.06211212667655827135426704899146*x1^5-0.052163923784865229293927768594585*x1^4+0.081140591279098853161322324467619*x1^3-0.077950826976361506370771792262531*x1^2-0.93405193717951635890006656381956*x1-0.020276638443425185065471794132463*x2^6-0.026242945812779666647784893029893*x2^5-0.12538335329298685993926198989357*x2^4+1.1492783526438503927113998770437*x2^3-0.11154140442300677915632434178406*x2^2-1.9297474994111820240094701262024*x2+0.81796822550108683191893987896037*x3^6+0.068966467071334303096108442332479*x3^5+1.102306576917776531621129265659*x3^4-4.4400325167154852390449804033778*x3^3-4.0498079986114728645585358890457*x3^2-0.33174966696044531910825270415444*x3-0.37205321834097124927831501395303];

gg = [1;1];
input = [0;gg*u1;gg*u2];
sym = [x2-x3^2; x3-x1^2+input(2); -x1-2*x2-x3+x2^3+input(3)];

%%
% V = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
% V = 1*x1^4+1*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; C0 = 5.862834287294065;
V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; C0 = 13.012408661959480;
%%
% f = [-x2-3/2*x1^2-1/2*x1^3; x1 - u];
% V = x1^2+x2^2+1*x1*x2;
% C0 = 3.00000001711;
% C0 = 2;
% V = 3*x1^2+2*x2^2+1*x1*x2;
% C0 = 2;
% C0 = 0.1;
%%
% f = [x2-x1
%     0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608
%     ];
% gg = 1;
% V = x1^2+x1*x2+x2^2;
% C0 = 100*8.109410321812767;
%%
dom = 8;
boundary_u = 100;
k = 4;
L_us = 4;
L_au = 4;
gamma = 2;
kk = 1;
i = 0; j = 0;
solh = C0 - V;
trace_Q = 0;
mm = 0;
TRACE = [];
Barrier = [];
%%
C1 = (x1-2)^2+(x2-1)^2+(x3-2)^2-1;
C2 = (x1+1)^2+(x2+2)^2+(x3+1)^2-1;
C3 = (x1+0)^2+(x2-0)^2+(x3-6)^2-9;
C4 = (x1+0)^2+(x2+0)^2+(x3+5)^2-9;
C5 = -x2+2;
C6 = (x1+0.5)^2+(x2+3.25)^2-1;
C7 = -x2 - 4;
C8 = -x1 + 1;
C = [C1;C2;C3;C4;C5;C6;C7;C8];
%%
figure(12);clf;hold on;view(-150, 30);
domain = [-dom dom -dom dom -dom dom];
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
%%%%%%%%%%%%%%%%%%%%%%%%
% while abs(double(trace_Q)-double(trace_Q1))>=epsi
% while 1
%
%     mm = mm+1; kk = 1;
%     fprintf('The whole Iteration time is:  %d\n  ',mm);

for i = 1:30
    fprintf('i=%6.0f\n',i);
    if kk == 0
        break
    else
        [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1_3D(f,k,solh,V,mm,gamma,gg,L_au,boundary_u);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if kk == 0
        break
    else
        [solh,trace_Q,kk]=sos_function_2_3D(f,k,SOLu1,SOLu2,SOL1,SOL2,gamma,mm,V,C,dom,gg,L_us);
        TRACE = [TRACE; double(trace_Q)];
        Barrier = [Barrier; solh];
    end
    %%%%%%%
end
xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);view(-150, 30);

%     solu =sos_function_3(k,solh,gamma,mm,V,C,V0);
%%
%     axis(domain)
%     kk = 1;
%     while kk == 1
%         j = j + 1;
%         fprintf('j=%6.0f\n',j);
%
%         record = solh;
%         record_Q = trace_Q;
%         [SOLu,SOL1,SOL2] = sos_function_4(k,solh,V,mm,gamma);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [solh, trace_Q, Q, kk]=sos_function_5(k,SOLu,SOL1,SOL2,gamma,mm,V,C);
%         if kk == 0
%             solh = record;
%             trace_Q = record_Q;
%         end
%     %%%%%%%
%     end

%     fprintf('The second round is end.====== \');
% %%
% end