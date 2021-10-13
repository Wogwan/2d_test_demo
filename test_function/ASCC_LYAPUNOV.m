pvar x1 x2;

plotVF_validate;hold on;
x = [x1;x2];
f = [x2; 
    (1-x1^2)*x2-x1];

deg = 2;
[Bs, gam] = barrier_test(deg);

V = gam - Bs;
% V = Bs - gam;
% V = Bs;

%%
dom = 4;
domain = [-dom dom -dom dom];

[~,~]=pcontour(V,double(C0),domain,'g'); hold on;             % Plot the original Lyapunov sublevel set
Vdot = jacobian(V,x)*f;
[~,~]=pcontour(Vdot,0,domain,'r'); hold on;


episC = 1e-8;                                                 % The condition to stop searching sub-level set C         
C0 = 0.01;                                                    % The initial value for searching sub-level set C
C = 0.1;                                                      % The initial value for searching sub-level set C                                                       % The initial value for recording steps for searching Barrier certificate
% domain = [-3 3 -3 3];
% domain = [-1.5 1.5 -1.5 1.5];                                 % The condition to stop searching barrier certificate
C_L_factor = 2;

%%
while abs(double(C)-double(C0))>=episC
    k = k + 1;
    fprintf('K =   %d\n ',k);
    if k==0
        L = search_LC_pvar_main(f,V,C0,C_L_factor);
    else
        C0 = C;
        L = search_LC_pvar_main(f,V,double(C),C_L_factor);
    end
    C = search_C_pvar_main(f,V,L)
end

plotVF_validate;
[~,~]=pcontour(V,double(C),domain,'b'); hold on;             % Plot the original Lyapunov sublevel set


%%
% C = double(C)
% h_t = C - V;                                                 % Set up the initial barrier certificate
% trace_CV = 0.1;
% trace_Q1 = 0;
% trace_Q = trace_CV;
% kkk = [];
% 
% %     while abs(double(trace_Q1)-double(trace_Q))>epsi
% while double(trace_Q)-double(trace_Q1)>=epsi
%     mm = mm+1;
%     fprintf('The Iteration time is:  %d\n  ',mm);
% 
%     if double(trace_Q)~=double(trace_CV)
%         trace_Q1 = trace_Q;
%     end
% 
%     %==========================================================================
%     % Algoritm step2: Check for the searching factor L1 and L2
%     [L1, L2] = search_L_pvar_main(f,V,h_t,L1_factor,L2_factor);
% %         [L1, L2] = search_L_pvar_main(f,V,h_t,L1_L2_factor,L1_L2_factor);
%     %==========================================================================
%     % Algoritm step3: Check for the qualified Barrier certificate
%     [h_t, trace_Q] = search_TracceQ_pvar_main(f,V,L1,L2,h_factor);
%     kkk = [kkk; double(trace_Q)];
%     % Conditionally plotting to record the barrier certificate
% 
%     if mod(mm,10)==1 || mod(mm,10)==2
%         [~,~]=pcontour(h_t,0,domain,'r'); hold on;          % Plot the original Lyapunov sublevel set
%     elseif mod(mm,10)==3 || mod(mm,10)==4
%         [~,~]=pcontour(h_t,0,domain,'m'); hold on;          % Plot the original Lyapunov sublevel set
%     elseif mod(mm,10)==5 || mod(mm,10)==6
%         [~,~]=pcontour(h_t,0,domain,'b'); hold on;          % Plot the original Lyapunov sublevel set
%     elseif mod(mm,10)==7 || mod(mm,10)==8
%         [~,~]=pcontour(h_t,0,domain,'c'); hold on;          % Plot the original Lyapunov sublevel set        
%     else
%         [~,~]=pcontour(h_t,0,domain,'g'); hold on;
%     end
%     % if mod(mm,20)==0
%     %     pause;
%     % end
%     refreshdata; drawnow; xlim([-1.2 1.2]); ylim([-1.2 1.2]); hold on;    
% end
