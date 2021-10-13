hold on;
syms x1 x2;
format long

%==========================================================================
% Initialize vector field Variable from %dotx2%
pvar x1 x2;
x = [x1;x2];
f = [x2-x1;
0.025359729472797157960711306552511*x1^6-0.01217009723261110433735733817297*x1^5-0.28561881441095026525735902480476*x1^4+0.0053102265416614026545590095169987*x1^3+x1^2*x2+0.17790534876627532116325861958709*x1^2-0.47421424943149480720390916606751*x1+0.083945383319193114801670674296474*x2^4+0.1356027624685912646995689101459*x2^3+0.071481552423721034239534333210031*x2^2-0.054409234652394423970012127256268*x2-0.000099625659358892892925041451235302
];
% f = [sys2pvar(dotx2(1));
%     sys2pvar(dotx2(2))];                                      % Convert the doxt3 from syms form to pvar form for SOS analysic
%     V = 6*x1^2+6*x2^2-8*x1*x2;                                  % 2 DEGREE Existed sublevel C and BC
%     V = -2*x1^2*x2^2+2*x1^4+6*x2^4+6*x1^2+6*x2^2-8*x1*x2;       % 4 DEGREE Existed sublevel C with 2 DEGREE BC
%     V = -2*x1^2*x2^2+2*x1^4+6*x2^4+6*x1^2+6*x2^2-4*x1*x2;       % 4 DEGREE Existed sublevel C with 4 DEGREE BC
%     V = 1*x1^2*x2^2+1*x1^4+6*x2^4+6*x1^2+6*x2^2-4*x1*x2;        % 4 DEGREE Existed sublevel C with 4 DEGREE BC
%     V = 6*(x1^4+2*x2^4-x1^2*x2^2)+2*(x1^2+x2^2-x1*x2);            % 4 DEGREE Existed sublevel C with 4 DEGREE BC
% V = 1*(2*x1^4+0.5*x2^4-x1^2*x2^2)+1*(1*x1^2+1*x2^2-0.5*x1*x2);
V = 1.279178112826987e-10*x1^4-4.298002103328788e-11*x1^2*x2^2+9.212657669110298e-13*x2^4+2.391619639114033e-12*x1^2+2.827357775433096e-12*x2^2;

episC = 1e-4;                                                 % The condition to stop searching sub-level set C         
epsi = 1e-5;                                                  % The condition to stop searching barrier certificate
C0 = 1e-2;                                                    % The initial value for searching sub-level set C
C = 0.1;                                                      % The initial value for searching sub-level set C
k = 0;                                                        % The initial value for recording steps for sub-level set C
km = 0;
mm = 0;                                                       % The initial value for recording steps for searching Barrier certificate
%     domain = [-3 3 -3 3];
domain = [-1.5 1.5 -1.5 1.5];                                 % The condition to stop searching barrier certificate
C_L_factor = 6;
L1_factor = 4;
L2_factor = 4;
h_factor = 6;

[~,~]=pcontour(V,double(C0),domain,'c'); hold on;             % Plot the original Lyapunov sublevel set
Vdot = jacobian(V,x)*f;
[~,~]=pcontour(Vdot,0,domain,'k'); hold on;


%==========================================================================

%     Algoritm step1: Check for the Sublevel set C
while double(C)-double(C0)>=episC
    k = k + 1;
    fprintf('K =   %d\n ',k);
    if k==0
        L = search_LC_pvar_main(f,V,C0,C_L_factor);
    else
        C0 = C;
        [L, k_mark] = search_LC_pvar_main(f,V,double(C),C_L_factor);
    end
    if k_mark == 1 && k == 1
        break
    end
    C = search_C_pvar_main(f,V,L);
    refreshdata; drawnow; 
end


[~,~]=pcontour(V,double(C),domain,'b'); hold on;             % Plot the original Lyapunov sublevel set

%     xlim([-2 2]);ylim([-2 2]);hold on;
xlim([-1.5 1.5]); ylim([-1.5 1.5]); hold on; 

% Compute the Original Barrier certificate
C = double(C)
h_t = C - V;                                                 % Set up the initial barrier certificate
trace_CV = 0.1;
trace_Q1 = 0;
trace_Q = trace_CV;
kkk = [];

%     while abs(double(trace_Q1)-double(trace_Q))>epsi
while double(trace_Q)-double(trace_Q1)>=epsi
    mm = mm+1;
    fprintf('The BC Iteration time is:  %d\n  ',mm);

    if double(trace_Q)~=double(trace_CV)
        trace_Q1 = trace_Q;
    end

    %==========================================================================
    % Algoritm step2: Check for the searching factor L1 and L2
    [L1, L2] = search_L_pvar_main(f,V,h_t,L1_factor,L2_factor);
%         [L1, L2] = search_L_pvar_main(f,V,h_t,L1_L2_factor,L1_L2_factor);
    %==========================================================================
    % Algoritm step3: Check for the qualified Barrier certificate
    [h_t, trace_Q] = search_TracceQ_pvar_main(f,V,L1,L2,h_factor);
    kkk = [kkk; double(trace_Q)];
    safe_BC_pvar = [safe_BC_pvar; h_t];
    % Conditionally plotting to record the barrier certificate

    if mod(mm,10)==1 || mod(mm,10)==2
        [~,~]=pcontour(h_t,0,domain,'r'); hold on;          % Plot the original Lyapunov sublevel set
    elseif mod(mm,10)==3 || mod(mm,10)==4
        [~,~]=pcontour(h_t,0,domain,'m'); hold on;          % Plot the original Lyapunov sublevel set
    elseif mod(mm,10)==5 || mod(mm,10)==6
        [~,~]=pcontour(h_t,0,domain,'b'); hold on;          % Plot the original Lyapunov sublevel set
    elseif mod(mm,10)==7 || mod(mm,10)==8
        [~,~]=pcontour(h_t,0,domain,'c'); hold on;          % Plot the original Lyapunov sublevel set        
    else
        [~,~]=pcontour(h_t,0,domain,'g'); hold on;
    end
    % if mod(mm,20)==0
    %     pause;
    % end
    refreshdata; drawnow; xlim([-1.2 1.2]); ylim([-1.2 1.2]); hold on;    
end