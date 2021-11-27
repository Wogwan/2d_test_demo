clear;tic;
pvar x1 x2 u htol epsi;
x = [x1;x2];
dom = 15;
%%
V = 0.04375278098272024*x1^4+0.002929993727326272*x1^3*x2-0.009153369646952632*x1^2*x2^2+0.1117128919832005*x1*x2^3+0.1337597502643151*x2^4-0.02068552493110127*x1^3-0.005144335364106218*x1^2*x2-0.06593999449113327*x1*x2^2-0.05150085197129949*x2^3+1.270472041414144*x1^2-0.04520082154572778*x1*x2+2.036327882785228*x2^2+0.0001597265387810012*x1+0.0001773175441663645*x2+18.32229824450674;
C0 = 41.61381252856700;
%%
k_u = 4; k_h = 4; L_us = 6; L_au = 4; gamma = 0;
k = ['r','g','b','m','c','k','y'];
%%
sol_B = C0 - V; solh = sol_B;
f = [x2-x1
    -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
    ];
gg = [1;1];
C1 = (x1+4)^2+(x2-6)^2-4;
C2 = (x1+3)^2+(x2+4)^2-4;
C3 = (x1-6)^2+(x2-0)^2-5;
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
B_controller_1 = 9.5200727248137635427838176838122*x1^2*x2^2-0.0019067415715053864094796765016326*x2-0.0027124382181331805306834237256908*x1-16.59753686623925972298820852302*x1*x2+63.590614823783475628715677885339*x1*x2^2+27.736623424839930152074884972535*x1^2*x2-57.175308245138516838323994306847*x1*x2^3-8.9090507766416919821494957432151*x1^3*x2-12.500188005127661483584233792499*x1^2+6.8887349167723126441842396161519*x1^3-6.1066090718013157356836018152535*x2^2-25.192936328667038026196678401902*x1^4+32.950912792629338810002082027495*x2^3-66.857596518457086176567827351391*x2^4+4858.2399374571323278360068798065;
[~,h31]=pcontour(B_controller_1,0,domain,'m'); hold on; h31.LineStyle = '-'; h31.LineWidth = 1.4;
axis(domain); TRACE = []; Barrier = []; Control = [];
%%
for j = 1:56
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