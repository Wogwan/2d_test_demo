clear;tic;
pvar x1 x2 u htol epsi;
x = [x1;x2];
%%
f = [x2-x1
    0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608
    ];
gg = [1;1];
dom = 5;
%%
V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; 
C0 = 59.903496442963558;
solh = C0 - V;
%%
k_u = 2;
k_h = 2;
L_us = 4;
L_au = 4;
% gamma = 0.1;
gamma = 2;
kk = 1;
i = 0; j = 0;
%%
C1 = (x1+3)^2+(x2-4)^2-1;
C2 = (x1+1)^2+(x2+4)^2-1;
C3 = (x1-3)^2+(x2-2)^2-1;
%%
C4 = (x1-2)^2+(x2-6)^2-1;
C5 = -x2+2;
C6 = (x1+0.5)^2+(x2+3.25)^2-1;
C7 = -x2 - 4;
C8 = -x1 + 1;
C = [C1;C2;C3;C4;C5;C6;C7;C8];
epsi = 1e-10;
trace_Q1 = 1;
trace_Q = 0;
mm = 0;
%%
figure(12);clf;hold on;
domain = [-dom dom -dom dom];
xlim([-dom dom]); ylim([-dom dom]); hold on;
[~,~]=pcontour(V,C0,domain,'b'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(C1,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(C2,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(C3,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(C4,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(C5,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(C6,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
% line([-dom,dom],[2,2],'linestyle','--');
axis(domain);
TRACE = [];
Barrier = [];
Control = [];
%%
% while abs(double(trace_Q)-double(trace_Q1))>=epsi
% while 1
%
%     mm = mm+1; kk = 1;
%     fprintf('The whole Iteration time is:  %d\n  ',mm);
% for i = 1:40
while 1
    i = i+1
    record_Q = trace_Q
    %%
    [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1(f,k_u,solh,V,gamma,gg,L_au);
    if kk == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2]];
    %%
    [solh,trace_Q,kk]=sos_function_2(i,f,k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us);
    if kk == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
end
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
toc