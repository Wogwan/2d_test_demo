clear;tic;
pvar x1 x2 u htol epsi;
x = [x1;x2];
dom = 15;
%%
V = 2.3527456760855520734310175612336*x1 + 2.0766624620968969772150103381136*x2 - 0.073876898960538034399903040139179*x1^2*x2^2 + 0.85470114203918201578602520385175*x1*x2 - 0.25680547023473881962374321119569*x1*x2^2 - 0.18352075410473428496160863687692*x1^2*x2 + 0.099754537825916755888755460546236*x1*x2^3 + 0.062942193715316535618242710370396*x1^3*x2 + 2.2453256315204943582841679017292*x1^2 + 0.35402259565560911802606369747082*x1^3 + 1.8873274087590221625987396691926*x2^2 + 0.26455550851938752776604246719216*x1^4 + 0.33022680287310124391808585642138*x2^3 + 0.2737765847508267236243284514785*x2^4 + 24.805093523509299302531871944666;
C0 = 48.09829077303992;
%%
k_u = 4; k_h = 4; L_us = 4; L_au = 6; gamma = 0;
k = ['r','g','b','m','c','k','y'];
%%
sol_B = C0 - V; solh = sol_B;
f = [x2-x1
    -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
    ];
gg = [1;1];
C1 = (x1+4)^2+(x2-5)^2-4;
C2 = (x1+0)^2+(x2+5)^2-4;
C3 = (x1-5)^2+(x2-0)^2-4;
C = [C1;C2;C3];
trace_Q1 = 1; trace_Q = 0; mm = 0; kk = 1; i = 0;
%%
figure_id = 12; figure(figure_id);clf;hold on;
domain = [-dom dom -dom dom]; xlim([-dom dom]); ylim([-dom dom]); hold on;
[~,~]=pcontour(V,C0,domain,'b'); hold on;
[~,~]=pcontour(C1,0,domain,'k'); hold on;
[~,~]=pcontour(C2,0,domain,'k'); hold on;
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'r');            % Plot the original Lyapunov sublevel set
end
B_1 = 1.3576759131504796407341473241104*x1*x2^2 - 3.0578816793752339187051347835222*x2 - 1.0488854530197364578469887419487*x1^2*x2^2 - 1.5393659377436941237959899808629*x1*x2 - 2.7790279263382875463150867290096*x1 + 0.99353551144585727783464790263679*x1^2*x2 - 1.0487971515579199710543889523251*x1*x2^3 - 1.0487550290522178464414082554867*x1^3*x2 - 0.59598573150985978408300525188679*x1^2 + 0.12189338665773290226734815178133*x1^3 - 0.4054861188650786818499227592838*x2^2 - 1.0495083410765784215357143693836*x1^4 + 0.77761775628344664834656896346132*x2^3 - 1.0494081234319465600890453060856*x2^4 + 89.818684691415668908121006097645;
[~,h31]=pcontour(B_1,0,domain,'m'); hold on; h31.LineStyle = '-'; h31.LineWidth = 1.4;
axis(domain); TRACE = []; Barrier = []; Control = [];
%%
for j = 1:80
% while 1
    i = i+1
    record_Q = trace_Q
    %%
    [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1(f,k_u,L_au,solh,V,gamma,gg);
    if kk == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2]];
    %%
    [solh,trace_Q,kk] = sos_function_2(i,f,k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id);
    if kk == 0
        break
    end
    figure(figure_id+1);clf;hold on;
    [~,h31]=pcontour(B_1,0,domain,'m'); hold on; h31.LineStyle = '-'; h31.LineWidth = 1.4;
    [~,~]=pcontour(V,C0,domain,'r'); hold on;
    [~,~]=pcontour(C1,0,domain,'k'); hold on;
    [~,~]=pcontour(C2,0,domain,'k'); hold on;
    [~,~]=pcontour(C3,0,domain,'k');
    if mod(i,7) == 0
        [~,~]=pcontour(solh,0,domain,k(7)); hold on;             % Plot the original Lyapunov sublevel set
    else
        [~,~]=pcontour(solh,0,domain,k(mod(i,7))); hold on;             % Plot the original Lyapunov sublevel set
    end
    refreshdata; drawnow;
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
end
fprintf('Optimal B(x) is \n%s \n\n',char(vpa(p2s(Barrier(end)))));
