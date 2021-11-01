% function [V,kk]=sos_optimal_v(f,gg,k,B,u,C)
clear;
pvar x1 x2 Vtol1 Vtol2;
x = [x1;x2];
dom = 50; domain = [-dom dom -dom dom];
%%
C1 = (x1+3)^2+(x2-4)^2-1;
C2 = (x1+3)^2+(x2+4)^2-1;
C3 = (x1-3)^2+(x2-2)^2-1;
C = [C1;C2;C3];
V0 = x1^2+x1*x2+x2^2;
C0_0 = 8.231821142289510;
B = -0.9580335390329583*x1^2+8.882107481271564e-05*x1*x2-3.30960108199112*x2^2+8.362588894146379e-05*x1+0.0002438801466290602*x2+37.66550908125493;
figure(14);clf;hold on;
% [~,~]=pcontour(V0,C0,domain,'y'); 
% [~,~]=pcontour(B,0,domain,'b');  
% [~,~]=pcontour(C(1),0,domain,'r');  
% [~,~]=pcontour(C(2),0,domain,'r');  
% [~,~]=pcontour(C(3),0,domain,'r');  
% [a1,b1] = coeffs(p2s(B));
% C0 = vpa(a1(end));
% b1(end)=[];
% a1(end)=[];
% B0 = s2p(a1*b1');
%%
f = [x2; -x1-x2*(1-x1^2)];
gg = [1;1];
u1 = -0.263205441705978*x1^2+1.861977417622602*x1*x2-1.956387277005004*x2^2-4207.835999213985*x1+5125.655463006833*x2-0.04136749982010112;
u2 = 0.3898848123891965*x1^2-1.961131454038698*x1*x2+2.010549285003059*x2^2+2951.536581729438*x1-4080.561465971713*x2+0.03995820559230309;
sys = [f(1)+gg(1)*u1; f(2)+gg(2)*u2];
kk = 1; gamma = 0;
V_degree = 2; l_au = 2;
l_us = 2; l_input = 2;
%%
[V,vc] = polydecvar('v_w',monomials(x,0:V_degree));
% [V,vc] = polydecvar('v_w',monomials(x,[1 V_degree]));
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:l_au));
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:l_au)); 
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:l_us)); 
[L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:l_us)); 
[L5,L5_Q] = sosdecvar('L5_w',monomials(x,0:l_us)); 
% [u1,u1_Q] = polydecvar('u1_w',monomials(x,0:l_input));
% [u2,u2_Q] = polydecvar('u2_w',monomials(x,0:l_input));
%%
hdot = jacobian(B, x1)*(f(1)+gg(1)*u1)+ jacobian(B, x2)*(f(2)+gg(2)*u2);
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
%% Constraint 1
% pcr_1 = L1 >= 0;
% pcr_2 = L2 >= 0;
% pcr_6 = L3 >= 0;
% pcr_7 = L4 >= 0;
% pcr_8 = L5 >= 0;
% pcr_31 = -V+C(1)*L3 >= 0;
% pcr_32 = -V+C(2)*L4 >= 0;
% pcr_33 = -V+C(3)*L5 >= 0;
% pconstr_1 = V-L1*B >= 0;
% pconstr_2 = Vdot+L2*B-Vtol >= 0;
% pconstr_3 = Vtol >= 0;
% % pconstr = [pcr_1; pcr_2; pcr_6; pcr_7; pcr_8; pcr_31; pcr_32; pcr_33; pconstr_1; pconstr_2; pconstr_3];
% % pconstr = [pcr_1; pcr_2; pcr_6; pcr_7; pcr_8; pconstr_1; pconstr_2; pconstr_3;pcr_31; pcr_32; pcr_33];
% % pconstr = [pconstr_1; pconstr_2; pconstr_3];
% pconstr = [pcr_1;pcr_2;pcr_6;pconstr_1;pconstr_2;pconstr_3;pcr_31];
%% Constraint
pcr_1 = L1 >= 0;
pcr_2 = L2 >= 0;
pcr_3 = L3 >= 0;
pcr_4 = V-C(1)*L3 >= 0;
pconstr_1 = V-L1*B >= 0;
% pconstr_1 = V-L1*B-Vtol1 >= 0;
pconstr_2 = -Vdot-L2*B+gamma*B-Vtol2 >= 0;
% pconstr_2 = Vdot+L2*B+gamma*B-Vtol2 >= 0;
% pconstr_2 = Vdot+L2*B+gamma*B >= 0;
% pconstr_3 = Vtol1 >= 0;
pconstr_4 = Vtol2 >= 0;
% pconstr = [pcr_1;pcr_2;pcr_6;pconstr_1;pconstr_2;pcr_31];
% pconstr = [pcr_1;pcr_2;pcr_3;pcr_4;pconstr_1;pconstr_2;pconstr_3;pconstr_4];
pconstr = [pcr_1;pcr_2;pcr_3;pcr_4;pconstr_1;pconstr_2;pconstr_4];
%% Set objection
obj = Vtol2;
% obj = Vtol1+Vtol2;
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,obj,opts);
% [info,dopt] = sosopt(pconstr,x,opts);

% Create output
if info.feas
    V = subs(V,dopt)
    
    % 16.6729
else
    kk = 0;
    V  = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end

%%
xlim([-dom dom]);ylim([-dom dom]);
% V = -0.0002688249895716166*x1^2+4.782656596144253e-05*x1*x2+0.001035911584956916*x2^2+6.53143778574437e-05*x1+9.265590495178159e-05*x2+722.2262195419746;
figure(15);clf;hold on;
cc = [1,2,3];
for i = 0:16
    k = ['r','g','b','m','c','k','y'];
    if mod(i,7) == 0
        [~,~]=pcontour(V,cc(1)*i,domain,k(7)); hold on;
        [~,~]=pcontour(V,cc(2)*i,domain,k(7)); hold on;
        [~,~]=pcontour(V,cc(3)*i,domain,k(7)); hold on;             
    else
        [~,~]=pcontour(V,cc(1)*i,domain,k(mod(i,7))); hold on;      
        [~,~]=pcontour(V,cc(2)*i,domain,k(mod(i,7))); hold on;      
        [~,~]=pcontour(V,cc(3)*i,domain,k(mod(i,7))); hold on;      
    end
end
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
% figure(16);clf;dom = 100; domain = [-dom dom -dom dom];
for i = 0:16
    k = ['r','g','b','m','c','k','y'];
    if mod(i,7) == 0
        [~,~]=pcontour(Vdot,0,domain,k(7)); hold on;             % Plot the original Lyapunov sublevel set
    else
        [~,~]=pcontour(Vdot,0,domain,k(mod(i,7))); hold on;             % Plot the original Lyapunov sublevel set
    end
end
% end
%% TEST FOR Sublevel Set
[a1,b1] = coeffs(p2s(V));
ccc = double(vpa(a1(end)));
C0 = ccc;
cc = ccc+0.1;
k = 2;
k_l = 4;
kk = 1;
%%
figure(14);hold on;
dom = 8; domain = [-dom dom -dom dom];
xlim([-dom dom]);ylim([-dom dom]);
[~,~]=pcontour(V0,C0_0,domain,'y'); 
[~,~]=pcontour(B,0,domain,'b'); 
[~,~]=pcontour(C(1),0,domain,'r');           
[~,~]=pcontour(C(2),0,domain,'r');                     
[~,~]=pcontour(V,ccc,domain,'g');      
solU = [];
v_c = [];
iter = 1;

%%
while abs(double(cc)-double(C0)) >= 1e-5
    iter = iter + 1;
    if iter ~= 1
        C0 = cc;
    end
    [solu1,solu2,solL,kk]= sos_optimal_v2(f,gg,k,k_l,V,C0);
    if kk == 0
        break
    end
    [cc,kk,solu1,solu2] = sos_optimal_v3(f,gg,k,k_l,V,C,dom,solL,ccc);
    v_c = [v_c; double(cc)]
    solU = [solU;[solu1,solu2]];
    if kk == 0
        figure(11);hold on;
        [~,~]=pcontour(V,v_c(end),domain,'g');  
        break
    end
end
toc