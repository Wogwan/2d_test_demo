clear;
load BC_test;
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
V_Sel = [];
C0_Sel = [];
for i = 1:length(A)
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
    u1 = A(i,1);
    u2 = A(i,2);
    B = A(i,3);
    %%
    c0 = 1; cc = 1.1; epsi = 1e-9; figure_id = 111; dom = 10; domain = [-dom dom -dom dom];
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
    %%
    [V, kk] = sos_optimal_V1_update(f,gg,B,u1,u2,V_au,V_us,V_degree,C,gamma);
    %%
    if kk == 0
        fprintf('Suitable Lyapunov function can not find.======\n');
    end
    
    %%
    figure(figure_id);hold on;
    [~,~]=pcontour(V,0,domain,'r');
   
end

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