clear;tic;
pvar x1 x2 u htol epsi;
x = [x1;x2];
dom = 15;
%%
f = [x2-x1
    -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
    ];
gg = [1;1];
%%
% V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
% C0 = 1.0e+02*1.237395136812358;  % Threshold 1e-9
% C0 = 1.0e+02*1.237395136880424; % Threshold 1e-6
% input = [gg(1)*u1;gg(2)*u2];
%%
V = 1*x1^4+1*x2^4+1*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
C0 = 1.0e+02*1.714593951979031;
%%
C1 = (x1+5)^2+(x2-8)^2-6;
C2 = (x1+7)^2+(x2+7)^2-6;
C3 = (x1-6)^2+(x2-0)^2-6;
C = [C1;C2;C3];
trace_Q1 = 1; trace_Q = 0;
mm = 0; kk = 1; i = 0;
%%
% f = [x2-x1;
%     0.46382156330341729589940137480476*x1^4-0.21219821363526560116570580989732*x1^3+0.35000677410250802565647474138031*x1^2+x1*x2-0.27258976376810471290805063896793*x1+0.15407267311901634565529661813343*x2^4+0.94268273772394006737584959410015*x2^3+4.1277331008261466394060335005634*x2^2-1.7098389065959235244562819389103*x2+0.011510115809394316777058975276304
%     ];
% gg = [1;1];
% V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
% if length(C) == 3
%     C0 = 123.7395131302175;
% else
%     C0 = 0.901712926410675; % V degree 2
% end
% V = x1^2+x1*x2+x2^2;
% C0 = 12.800271326769639; % V degree 2
%%
sol_B = C0 - V;
solh = sol_B;
%%
% k_u = 4; k_h = 4;
% L_us = 4; L_au = 2; gamma = 0;
%% Test
k_u = 4; k_h = 4; L_us = 4; L_au = 4;
gamma = 0;
%%
figure_id = 12;
figure(figure_id);clf;hold on;
domain = [-dom dom -dom dom];
xlim([-dom dom]); ylim([-dom dom]); hold on;
[~,~]=pcontour(V,C0,domain,'b'); hold on;
[~,~]=pcontour(C1,0,domain,'k'); hold on;
[~,~]=pcontour(C2,0,domain,'k'); hold on;
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'k');            % Plot the original Lyapunov sublevel set
end
axis(domain); TRACE = [];
Barrier = []; Control = [];
%%
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
    %     [solh,trace_Q,kk]=sos_function_2(i,f,k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us,sol_B);
    [solh,trace_Q,kk] = sos_function_2(i,f,k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id);
    if kk == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
end

A = [];
for i = 1:length(Barrier)
    A = [A; [Control(i,:) Barrier(i)]];
end