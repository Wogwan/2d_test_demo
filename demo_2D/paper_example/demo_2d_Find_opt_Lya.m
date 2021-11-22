clear;
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
f = [x2-x1
    -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
    ];
gg = [1; 1];
sys = [f(1)+gg(1)*u1; f(2)+gg(2)*u2];
%%
C1 = (x1+4)^2+(x2-6)^2-4;
C2 = (x1+3)^2+(x2+4)^2-4;
C3 = (x1-6)^2+(x2-0)^2-5;
C = [C1;C2;C3];
V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
% C0 = 1.0e+02*1.237395136955980;
C0 = 1.0e+02*1.237395136880424;
%% Start to compute the control barrier function
%% Threshold 1e-9 Original version
% u1 = 6.100196267619156e-08*x1^4-1.546015265265402e-07*x1^3*x2+1.113095358488676e-06*x1^2*x2^2+1.092215740139554e-06*x1*x2^3-3.511916361791863e-07*x2^4-87.2633504212407*x1^3+70.4373421754175*x1^2*x2-450.7021749132224*x1*x2^2+97.44006698690981*x2^3-14.37926152094848*x1^2+125.4238863613426*x1*x2-223.3072193661567*x2^2-100.5122508004279*x1+22.81227168299663*x2-0.006327307808920016;
% u2 = 0.2360635815306517*x1^4+6.985184784813693e-08*x1^3*x2+8.621740491180531e-08*x1^2*x2^2-3.210134156589625e-07*x1*x2^3-0.03729669708706964*x2^4+167.5019986103727*x1^3-176.9870680058685*x1^2*x2+138.7018787583334*x1*x2^2-43.89829329803597*x2^3+33.77638723860747*x1^2-92.7653076474012*x1*x2+32.85413854628985*x2^2+60.28003840340249*x1-91.83277810273101*x2-0.2417279416649734;
% B = -47.45989365990555*x1^4-16.2666888328201*x1^3*x2+18.75887722919024*x1^2*x2^2-107.0338306842238*x1*x2^3-125.0019692594601*x2^4+10.85612639957994*x1^3+52.65090970523237*x1^2*x2+117.9496617350447*x1*x2^2+58.03043337506556*x2^3-21.51146693599407*x1^2-27.65196501886886*x1*x2-9.46303451193095*x2^2-0.006725081215624065*x1-0.004566790324176556*x2+8969.094719261169;
%% Threshold 1e-6 Latest version
u1 = 2.906220856916317e-08*x1^4+3.94038660700824e-08*x1^3*x2+4.191725481939133e-07*x1^2*x2^2+3.661557007086112e-07*x1*x2^3-4.243496200651684e-07*x2^4-57.65528482511168*x1^3+21.84790702043916*x1^2*x2-219.6672176479056*x1*x2^2+41.10348829759877*x2^3-13.79278726113727*x1^2+49.4041054600597*x1*x2-166.319805170988*x2^2-100.5440239810693*x1+17.10982649929812*x2-0.004597195361128125;
u2 = 0.2360640082093554*x1^4+9.85864162586925e-08*x1^3*x2+8.169529216323998e-08*x1^2*x2^2-2.111889208524325e-07*x1*x2^3-0.03729668002161401*x2^4+77.95805512315057*x1^3-88.55295447165489*x1^2*x2+68.32402608201978*x1*x2^2-28.1361460228466*x2^3+31.75121039639564*x1^2-47.70000074122429*x1*x2+19.59605552885746*x2^2+48.69554776307634*x1-75.76839093443586*x2-0.2392243380101325;
B = -25.19293632866704*x1^4-8.909050776641692*x1^3*x2+9.520072724813764*x1^2*x2^2-57.17530824513852*x1*x2^3-66.85759651845709*x2^4+6.888734916772313*x1^3+27.73662342483993*x1^2*x2+63.59061482378348*x1*x2^2+32.95091279262934*x2^3-12.50018800512766*x1^2-16.59753686623926*x1*x2-6.106609071801316*x2^2-0.002712438218133181*x1-0.001906741571505386*x2+4858.239937457132;
%%
c0 = 1; cc = 1.1; epsi = 1e-9; figure_id = 111; dom = 20; domain = [-dom dom -dom dom];
%%
figure(figure_id);clf;hold on;
[~,~]=pcontour(C(1),0,domain,'r');
[~,~]=pcontour(C(2),0,domain,'r');
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'r');
elseif length(C) == 4
    [~,~]=pcontour(C(3),0,domain,'r');
    [~,~]=pcontour(C(4),0,domain,'r');
else
    fprintf('The constraint number does not match.======\n');
end
hold on; [~,~]=pcontour(B,0,domain,'g');
[~,~]=pcontour(V0,C0,domain,'m');
%% Hyperparameters of the SOSP @ CBF -> V
V_us = 4; V_au = 4; V_degree = 4; gamma = 0; k_u_V = 4; k_l_au = 4;
% V_us = 6; V_au = 6; V_degree = 6; gamma = 0; k_u_V = 6; k_l_au = 6;
% V_us = 8; V_au = 8; V_degree = 8; gamma = 0; k_u_V = 8; k_l_au = 8;
%%
kk = 1; OO = 0;
dom_2 = 10000; domain_2 = [-dom_2 dom_2 -dom_2 dom_2]; N_Lya = [];
%%
[V, kk] = sos_optimal_V1(f,gg,B,u1,u2,V_au,V_us,V_degree,C,gamma);
%%
if kk == 0
    fprintf('Suitable Lyapunov function can not find.======\n');
end
N_Lya = [N_Lya;V];
%% TEST FOR Sublevel Set
[a1,b1] = coeffs(p2s(V));
ccc = double(vpa(a1(end)));
C0 = ccc; cc = ccc+1;

%%
figure(figure_id);hold on;
[~,~]=pcontour(V,double(C0),domain,'r');
solU = []; v_c = []; iter = 0;
%%
% while double(cc)-double(C0) >= epsi
while double(cc)-double(C0) >= 1e-5
    iter = iter + 1;
    if iter ~= 0
        C0 = cc;
    end
    [solL,kk]= sos_optimal_v2(f,gg,k_u_V,k_l_au,V,cc);
    if kk == 0
        break
    end
    [cc,kk,solu1,solu2] = sos_optimal_v3(f,gg,k_u_V,k_l_au,V,C,dom,solL,ccc,figure_id);
    v_c = [v_c; double(cc)]
    solU = [solU;[solu1,solu2]];
    if kk == 0
        figure(figure_id);hold on;
        [~,~]=pcontour(V,v_c(end),domain,'b');
        break
    end
end

%% Start to compute the control barrier function
c_b = v_c(end);
sol_B = c_b - N_Lya(end);
hold on;[~,~]=pcontour(sol_B,0,domain,'k');
V = N_Lya(end);

%%
%     i_count_1 = min(length(Control),length(Barrier));
%     if i_count_1 ~= 0
%         u1 = Control(i_count_1,1);
%         u2 = Control(i_count_1,2);
%         B = Barrier(i_count_1);
%     else
%         if isempty(Control)
%             u1 = 1;
%             u2 = 1;
%         else
%             u1 = Control(i_count_1,1);
%             u2 = Control(i_count_1,2);
%         end
%         B = 0;
%         fprintf('Barrier function can not find.======\n');
%     end

%     figure(figure_id+2);clf;
%     hold on;
%     cc = [1,2,3];
%     for i = 0:16
%         k_u_V = ['r','g','b','m','c','k','y'];
%         if mod(i,7) == 0
%             [~,~]=pcontour(V,cc(1)*i,domain_2,k_u_V(7)); hold on;
%             [~,~]=pcontour(V,cc(2)*i,domain_2,k_u_V(7)); hold on;
%             [~,~]=pcontour(V,cc(3)*i,domain_2,k_u_V(7)); hold on;
%         else
%             [~,~]=pcontour(V,cc(1)*i,domain_2,k_u_V(mod(i,7))); hold on;
%             [~,~]=pcontour(V,cc(2)*i,domain_2,k_u_V(mod(i,7))); hold on;
%             [~,~]=pcontour(V,cc(3)*i,domain_2,k_u_V(mod(i,7))); hold on;
%         end
%     end
%     Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
%     for i = 0:16
%         k_u_V = ['r','g','b','m','c','k','y'];
%         if mod(i,7) == 0
%             [~,~]=pcontour(Vdot,0,domain_2,k_u_V(7)); hold on;             % Plot the original Lyapunov sublevel set
%         else
%             [~,~]=pcontour(Vdot,0,domain_2,k_u_V(mod(i,7))); hold on;             % Plot the original Lyapunov sublevel set
%         end
%     end