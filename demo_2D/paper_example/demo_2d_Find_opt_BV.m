clear;tic;
pvar x1 x2 u htol epsi;
x = [x1;x2];
dom = 10;

%% Optimize without vtol2
% V = 0.02930787526537085*x1^4+0.0140334521199723*x1^3*x2-0.003427601802693766*x1^2*x2^2+0.09183813309170435*x1*x2^3+0.09869600181719627*x2^4-0.01153095893518936*x1^3-0.01486925374081773*x1^2*x2-0.05665519025240803*x1*x2^2-0.0531466040980413*x2^3+1.429403099563379*x1^2-0.2299306443098457*x1*x2+2.524970194719252*x2^2+0.0002143489465478355*x1+0.0003197678206660723*x2+1.532702817990033;
% C0 = 25.765249810193705; 
% k_u = 4; k_h = 4; L_us = 4; L_au = 6; gamma = 0;

%%
V = 0.04791793554535747*x1^4+0.005751975855669967*x1^3*x2-0.01465931654811854*x1^2*x2^2+0.1264746534640889*x1*x2^3+0.1435738234375542*x2^4-0.02475700871808066*x1^3-0.01677068421056823*x1^2*x2-0.07592877484460123*x1*x2^2-0.06910059119814008*x2^3+1.002549245831866*x1^2-0.02331031823434397*x1*x2+1.530088835906439*x2^2+0.0001997984676651864*x1+0.0002674429123316613*x2+7.669980250359338;
C0 = 28.430098121975416; 
k_u = 4; k_h = 4; L_us = 6; L_au = 4; gamma = 0;

%%
sol_B = C0 - V; solh = sol_B;
f = [x2-x1
    -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
    ];
gg = [1;1];
C1 = (x1+4)^2+(x2-6)^2-4;
C2 = (x1+3)^2+(x2+4)^2-4;
C3 = (x1-6)^2+(x2-0)^2-5;
C = [C1;C2;C3]; % C = [C1;C2];
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
B_controller_1 = -25.5286775994005*x1^4-8.826737822184512*x1^3*x2+9.97662491998584*x1^2*x2^2-57.03373052165595*x1*x2^3-66.6347601541668*x2^4+6.719666264862423*x1^3+27.82600314594399*x1^2*x2+63.24796692083677*x1*x2^2+32.27190868643367*x2^3-12.33737345606981*x1^2-16.10987157050511*x1*x2-5.888600774219452*x2^2-0.003646379657348746*x1-0.002562037970909413*x2+4830.062776216545;
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
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
end