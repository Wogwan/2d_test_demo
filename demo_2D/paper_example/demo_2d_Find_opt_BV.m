clear;tic;
pvar x1 x2 u htol epsi;
x = [x1;x2];
dom = 15;
%%
V = 0.6536260414188650447187001191196*x1 + 0.70808287458676077985586516660987*x2 + 0.11335671197243221697270598724572*x1^2*x2^2 - 0.36180728954161472943340527308465*x1*x2 - 0.028558208518183213209251647413112*x1*x2^2 - 0.042161977404113662459828049122734*x1^2*x2 + 0.13615092192557509687134142950526*x1*x2^3 + 0.11075242434548854264519945900247*x1^3*x2 + 0.87368437299271362039831956280977*x1^2 + 0.11464565489933541131417626957045*x1^3 + 0.87548828852503746134061657357961*x2^2 + 0.1018082983738275160146002917827*x1^4 + 0.037909772744884030759582316250089*x2^3 + 0.12344554380472853860606363696206*x2^4 + 18.850890322167181523127510445192;
C0 = 33.38322620290549;
%%
k_u = 4; k_h = 4; L_us = 4; L_au = 2; gamma = 0;
k = ['r','g','b','m','c','k','y'];
%%
sol_B = C0 - V; solh = sol_B;
f = [x2-x1
    -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
    ];
gg = [1;1];
C1 = (x1+4)^2+(x2-6)^2-4;
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
B_1 = 8.4371824851844099413256117259152*x1*x2^2 - 12.108574171358178261925786500797*x2 - 4.7933188908649331239075763733126*x1^2*x2^2 - 1.8005404644777436296010364458198*x1*x2 - 7.6948518605698605909992693341337*x1 + 7.4772107108042220602328598033637*x1^2*x2 - 4.7934462359877381132378104666714*x1*x2^3 - 4.7932375160978821782009617891163*x1^3*x2 + 3.6524473614994743542183641693555*x1^2 + 2.8271808409055445565627451287583*x1^3 + 2.0036370641438416839719138806686*x2^2 - 4.7929211951298178462366195162758*x1^4 + 7.6872791727953684315366444934625*x2^3 - 4.7953514892347675058204004017171*x2^4 + 494.16343828145528505046968348324;
[~,h31]=pcontour(B_1,0,domain,'m'); hold on; h31.LineStyle = '-'; h31.LineWidth = 1.4;
axis(domain); TRACE = []; Barrier = []; Control = [];
%%
% for j = 1:60
while 1
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
