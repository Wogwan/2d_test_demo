clear;
% load BC_test_2D;
% load BC_test_2D2;
% load R_Change1_BC_test_2D;
load R_Change1_BC_test_2D2;
% load R_Change1_BC_test_2D3;
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
k = ['r','g','b','m','c','k','y'];
dom = 50; domain = [-dom dom -dom dom];
V_Sel = []; C0_Sel = []; N_Lya = [];
figure_id = 111;
%%
for i = 1:length(A(:,1))
    i
    f = [x2-x1
        -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
        ];
    gg = [1; 1];
    sys = [f(1)+gg(1)*u1; f(2)+gg(2)*u2];
    %%
    %     C1 = (x1+5)^2+(x2-8)^2-6;
    %     C2 = (x1+7)^2+(x2+7)^2-6;
    %     C3 = (x1-6)^2+(x2-0)^2-6;
    C1 = (x1+5)^2+(x2-8)^2-6;
    C2 = (x1+7)^2+(x2+7)^2-6;
    C3 = (x1-6)^2+(x2+2)^2-6;
    C = [C1;C2;C3];
    V0 = 1*x1^4+1*x2^4+1*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
    %     C0 = 1.0e+02*1.714593951979031;
    C0 = 1.0e+2*2.121924848371969;
    %% Start to compute the control barrier function
    u1 = A(i,1);
    u2 = A(i,2);
    %     B = A(i,3);
    B = A(i,4);
    %%
    c0 = 1; cc = 1.1; epsi = 1e-5; N_Lya = [];
    %%
    figure(figure_id);clf;hold on;
    [~,~]=pcontour(C(1),0,domain,'r');
    [~,~]=pcontour(C(2),0,domain,'r');
    if length(C) == 3
        [~,~]=pcontour(C(3),0,domain,'r');
    else
        fprintf('The constraint number does not match.======\n');
    end
    [~,~]=pcontour(V0,C0,domain,'m');
    [~,~]=pcontour(B,0,domain,'B');
    %% Hyperparameters of the SOSP @ CBF -> V
    V_us = 4;V_au = 6;V_degree = 4;k_u_V = 4;k_l_au = 4;gamma = 0;
    %     V_us = 4;V_au = 4;V_degree = 4;k_u_V = 4;k_l_au = 4;gamma = 0;
    % V_us = 2;V_au = 4;V_degree = 2;k_u_V = 2;k_l_au = 4;gamma = 0;
    %     V_us = 6; V_au = 6; V_degree = 6; gamma = 0; k_u_V = 6; k_l_au = 6;
    % V_us = 8; V_au = 8; V_degree = 8; gamma = 0; k_u_V = 8; k_l_au = 8;
    %%
    kk = 1; OO = 0;
    %%
    %     [V, kk] = sos_optimal_V1_2D_V(f,gg,B,u1,u2,V_au,V_us,V_degree,C,gamma);
    [V, kk] = sos_optimal_V1(f,gg,B,u1,u2,V_au,V_us,V_degree,C,gamma);
    %%
    if kk == 0
        fprintf('53 Suitable Lyapunov function can not find.======\n');
    end
    N_Lya = [N_Lya;V];
    %% TEST FOR Sublevel Set
    [a1,b1] = coeffs(p2s(V));
    ccc = double(vpa(a1(end)));
    C0 = ccc; cc = ccc+0.5;
    solU = []; v_c = []; iter = 0;
    %%
    figure(figure_id);hold on;
    [~,~]=pcontour(V,double(C0),domain,'r');
    %%
    %     while abs(double(cc)-double(C0)) >= epsi
    while double(cc)-double(C0) >= epsi
        iter = iter + 1;
        if iter ~= 0
            C0 = cc;
        end
        %         [solL,kk]= sos_optimal_v2(f,gg,k_u_V,k_l_au,V,cc);
        [solL,kk]= sos_optimal_v2_origin(f,gg,k_u_V,k_l_au,V,cc);
        if kk == 0
            break
        end
        %         [cc,kk,solu1,solu2] = sos_optimal_v3(f,gg,k_u_V,k_l_au,V,C,dom,solL,ccc,figure_id);
        [cc,kk,solu1,solu2] = sos_optimal_v3_origin(f,gg,k_u_V,k_l_au,V,C,dom,solL,ccc,figure_id);
        v_c = [v_c; double(cc)];
        solU = [solU;[solu1,solu2]];
        if kk == 0
            figure(figure_id);hold on;
            [~,~]=pcontour(V,v_c(end),domain,'b');
            break
        end
        %%
        %         if mod(iter,7) == 0
        %             [~,~]=pcontour(V,double(cc),domain,k(7)); hold on;             % Plot the original Lyapunov sublevel set
        %         else
        %             [~,~]=pcontour(V,double(cc),domain,k(mod(iter,7))); hold on;             % Plot the original Lyapunov sublevel set
        %         end
        refreshdata; drawnow;
    end
    %% Start to compute the control barrier function
    c_b = v_c(end);
    sol_B = c_b - N_Lya(end);
    hold on;[~,~]=pcontour(sol_B,0,domain,'k');
    V = N_Lya(end);
    V_Sel = [V_Sel; V];
    C0_Sel = [C0_Sel; c_b];
    %%
    figure(figure_id+1);clf;hold on;
    [~,~]=pcontour(C(1),0,domain,'k');
    [~,~]=pcontour(C(2),0,domain,'k');
    if length(C) == 3
        [~,~]=pcontour(C(3),0,domain,'k');
    else
        fprintf('The constraint number does not match.======\n');
    end
    %     pcontour(V,0,domain,'B');
    pcontour(sol_B,0,domain,'g');
    xlim([-dom dom]); ylim([-dom dom]);
end
