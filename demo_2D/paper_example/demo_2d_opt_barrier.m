clear;tic;
pvar x1 x2 u htol epsi;
format long
x = [x1;x2];
dom = 10;
%%
f = [x2-x1
    -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
    ];
gg = [1;1];
%%
% V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; 
V = 1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1;
C0 = 57.23191153176421;
%%
C1 = (x1+4)^2+(x2-5)^2-4;
C2 = (x1+0)^2+(x2+5)^2-4;
C3 = (x1-5)^2+(x2-0)^2-4;
%%
C = [C1;C2;C3];
trace_Q1 = 1; trace_Q = 0; mm = 0; kk = 1; iter = 0;
k = ['r','g','b','m','c','k','y'];
%%
sol_B = C0 - V;
solh = sol_B;
%%
% k_u = 4; k_h = 4; L_us = 4; L_au = 2; gamma = 0;
k_u = 4; k_h = 4; L_us = 4; L_au = 2; gamma = 0;
%%
figure_id = 12;
figure(figure_id+1);clf;hold on;
figure(figure_id+2);clf;hold on;
figure(figure_id);clf;hold on;
domain = [-dom dom -dom dom];
xlim([-dom dom]); ylim([-dom dom]); hold on;
[~,~]=pcontour(V,C0,domain,'b'); hold on;
[~,~]=pcontour(C1,0,domain,'k'); hold on;
[~,~]=pcontour(C2,0,domain,'k'); hold on;
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'k');         
end
axis(domain); TRACE = [];
Barrier = []; Control = []; Barrier_plus = [];
%%
while 1
    iter = iter+1
    record_Q = trace_Q
    %%
    [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1(f,k_u,L_au,solh,V,gamma,gg);
    if kk == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2]];
    %%
    [solh,trace_Q,kk] = sos_function_2(iter,f,k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id);
    if kk == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
    %% Optimal the set
    kk = 1; OO = 0;
    %%
    if kk == 0
        fprintf('Advanced Barrier Function can not find.======\n');
    end
    figure(figure_id+1);clf;hold on;
    [~,~]=pcontour(V,C0,domain,'b'); hold on;
    [~,~]=pcontour(C1,0,domain,'k'); hold on;
    [~,~]=pcontour(C2,0,domain,'k'); hold on;
    [~,~]=pcontour(C3,0,domain,'k');
    [~,~]=pcontour(V,C0,domain,'r'); hold on;
    if mod(iter,7) == 0
        [~,~]=pcontour(solh,0,domain,k(7)); hold on;             % Plot the original Lyapunov sublevel set
    else
        [~,~]=pcontour(solh,0,domain,k(mod(iter,7))); hold on;             % Plot the original Lyapunov sublevel set
    end
    refreshdata; drawnow;
end
toc
A = [];
for iter = 1:length(Barrier)
    A = [A; [Control(iter,:) Barrier(iter)]];
end
%%
fprintf('Permissive B(x) is \n%s \n\n',char(vpa(p2s(A(end,3)))));
fprintf('Control Input u1(x) is \n%s \n\n',char(vpa(p2s(A(end,1)))));
fprintf('Control Input u2(x) is \n%s \n\n',char(vpa(p2s(A(end,2)))));