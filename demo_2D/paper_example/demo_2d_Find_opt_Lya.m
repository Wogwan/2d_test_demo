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
c0 = 1; cc = 1.1; epsi = 1e-6; epsi2 = 0.196;
figure_id = 111;
%%
dom = 15; domain = [-dom dom -dom dom];
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
[~,~]=pcontour(V0,c0,domain,'g');
%% Start to compute the control barrier function
c_b = 12.800276259538252;
sol_B = c_b - V0;
solh = sol_B;
%% Hyperparameters of the SOSP @ CBF
% b_k_u = 6; b_k_h = 4;
% L_us = 8; L_au = 8;
% gamma = 0;
b_k_u = 4; b_k_h = 4;
L_us = 8; L_au = 8;
gamma = 0;
kk = 1; j = 0;
TRACE = []; Barrier = []; Control = [];
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

u1 = -4.145115472236079e-08*x1^4+1.085346287651691e-07*x1^3*x2-1.063440665979916e-06*x1^2*x2^2-1.273284365624724e-06*x1*x2^3+3.356899144529559e-08*x2^4-68.10857648776478*x1^3+44.79257440458176*x1^2*x2-317.3794055004882*x1*x2^2+74.31186325734753*x2^3-13.24233704234717*x1^2+84.57833998544589*x1*x2-183.0342177990875*x2^2-96.11420842778074*x1+20.28668806951276*x2-0.005856492685127363;
u2 = 0.2360647594407567*x1^4+5.153703495909361e-08*x1^3*x2-4.439371098023446e-08*x1^2*x2^2+2.685979145088956e-07*x1*x2^3-0.0372967752260427*x2^4+116.9246410032892*x1^3-128.6976472804574*x1^2*x2+100.7124626933768*x1*x2^2-36.15376650850188*x2^3+29.89385165904772*x1^2-61.7623336896533*x1*x2+25.67702208750281*x2^2+53.96918394213625*x1-81.61150344848821*x2-0.2401018187221746;
B = -25.5286775994005*x1^4-8.826737822184512*x1^3*x2+9.97662491998584*x1^2*x2^2-57.03373052165595*x1*x2^3-66.6347601541668*x2^4+6.719666264862423*x1^3+27.82600314594399*x1^2*x2+63.24796692083677*x1*x2^2+32.27190868643367*x2^3-12.33737345606981*x1^2-16.10987157050511*x1*x2-5.888600774219452*x2^2-0.003646379657348746*x1-0.002562037970909413*x2+4830.062776216545;
hold on; [~,~]=pcontour(B,0,domain,'g');
C0 = 1.0e+02*1.237395136955980;
V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
[~,~]=pcontour(V0,C0,domain,'m');

%%

%% Hyperparameters of the SOSP @ CBF -> V
V_us = 4; V_au = 4; V_degree = 4;
gamma = 0;
kk = 1; OO = 0;
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
%%
dom_2 = 10000; domain_2 = [-dom_2 dom_2 -dom_2 dom_2]; N_Lya = [];
%%
[V, kk] = sos_optimal_V1(f,gg,B,u1,u2,V_au,V_us,V_degree,C,gamma);
%%
if kk == 0
    fprintf('Suitable Lyapunov function can not find.======\n');
end
N_Lya = [N_Lya;V];
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
%% TEST FOR Sublevel Set
[a1,b1] = coeffs(p2s(V));
ccc = double(vpa(a1(end)));
C0 = ccc; cc = ccc+0.1;
k_u_V = 4; k_l_au = 4; kk = 1;
%%
figure(figure_id);hold on;
[~,~]=pcontour(V,double(cc),domain,'c');
solU = []; v_c = []; iter = 0;
%%
while abs(double(cc)-double(C0)) >= epsi
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
