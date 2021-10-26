clear;

pvar x1 x2 u htol epsi;
x = [x1;x2];
%%%%%%%%%%%%%%%%%%%%%%%
%%
% f = [x2; -x1 + u];
% V = x1^2+x2^2+x1*x2;
%%
% f = [-x2-3/2*x1^2-1/2*x1^3; x1 - u];
% V = x1^2+x2^2+1*x1*x2;
% C0 = 3.00000001711;
% C0 = 2;
% V = 3*x1^2+2*x2^2+1*x1*x2;
% C0 = 2;
% C0 = 0.1;
%%
f = [0.1*x1^2+1*x2; 0.1*x1*x2-0.2*x1+(1+x1^2)*u];
V = 1*x1^2+1*x2^2+1*x1*x2;
C0 = 5.67547012045;

k = 2;
gamma = 1;
kk = 1;
i = 0; j = 0;
solh = C0 - V;
V0 = solh;
C1 = (x1-3)^2+(x2-1)^2-1;
C2 = (x1+3)^2+(x2+4)^2-1;
C3 = (x1+4)^2+(x2-5)^2-1;
C4 = (x1-2)^2+(x2-6)^2-1;
C5 = -x2+2;
C6 = (x1+0.5)^2+(x2+3.25)^2-1;
C7 = -x2 - 4;
C8 = -x1 + 1;
C = [C1;C2;C3;C4;C5;C6;C7;C8];

epsi = 1e-5;
trace_Q1 = 1;
trace_Q = 0;
mm = 0;
saved_BC = [];                                            % Set to store the learned Barrier Certificate function
saved_u = [];                                            % Set to store the learned Barrier Certificate function
solh_re = solh;

figure(11);clf;hold on;
dom = 25;
domain = [-dom dom -dom dom];
% domain = [-20 20 -20 20];
xlim([-dom dom]); ylim([-dom dom]); hold on;   
[~,~]=pcontour(solh,0,domain,'c'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(C1,0,domain,'r'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(C2,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(C3,0,domain,'b'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(C4,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(C5,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(C6,0,domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
% line([-dom,dom],[2,2],'linestyle','--');

axis(domain);

%%%%%%%%%%%%%%%%%%%%%%%%

% while abs(double(trace_Q)-double(trace_Q1))>=epsi
% while 1
%     
%     mm = mm+1; kk = 1;
%     fprintf('The whole Iteration time is:  %d\n  ',mm);

    while kk == 1
        i = i + 1;
        fprintf('i=%6.0f\n',i);
        record = solh;
        record_Q = trace_Q;
        [SOLu,SOL1,SOL2] = sos_function_1(k,solh,V,mm,gamma,dom);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [solh, trace_Q, Q, kk]=sos_function_2(f,k,SOLu,SOL1,SOL2,gamma,mm,V,C,dom);
        if kk == 0
            solh = record;
            trace_Q = record_Q;
        end
    %%%%%%%
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