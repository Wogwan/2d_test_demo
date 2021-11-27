% clc;
% load BC_test_Lya;
% load BC_test_Lya2;
% load test_1125_1_3D;
% load test_1125_3_3D;
% load figure_3D;
% clear;
KK = 6;
% load test_1126_3D_opt_BA_compute;
Barrier = Compute{KK};
%%
load test_1126_1_3D;
B1 = A(KK,4);
%%
pvar x1 x2 x3;
dom = 10; domain = [-dom dom -dom dom -dom dom];
C1 = (x1-3)^2+(x2-0)^2+(x3-0)^2-2;
C2 = (x1+3)^2+(x2+3)^2+(x3-3)^2-3;
C3 = (x1-0)^2+(x2-3)^2+(x3+0)^2-3;
C4 = (x1-3)^2+(x2-0)^2+(x3+3)^2-3;
% V = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
% C0 = 9.681159900631842;
% V = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
% C0 = 2.584683714740699;
k = ['r','g','m','c','k','y'];
%%
% for iter = 1:length(A)
%     B = A(iter,4);
%%
iter = 1;
Ba = Barrier{iter};
for iter = 1:length(Ba)
    B = Ba(iter,1);
    %%
    figure(iter);clf;hold on;
    phB= patch(pcontour3(B,0,domain,'B'));
    %     phV0= patch(pcontour3(V,double(C0),domain,'c')); set(phV0,'EdgeAlpha',1, 'FaceColor', 'none', 'EdgeColor', 'b' );
    phB1= patch(pcontour3(B1,0,domain,'B')); set(phB1,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',1);
    %%
    %     ph1= patch(pcontour3(C1,0,domain,'c')); set(ph1, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph2= patch(pcontour3(C2,0,domain,'c')); set(ph2, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph3= patch(pcontour3(C3,0,domain,'c')); set(ph3, 'FaceColor', 'none', 'EdgeColor', 'k' );
    ph4= patch(pcontour3(C4,0,domain,'c')); set(ph4, 'FaceColor', 'none', 'EdgeColor', 'k' );
    
    if mod(iter,6) == 0
        set(phB,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', k(6),'LineStyle','-','LineWidth',1);
    else
        set(phB,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', k(mod(iter,6)),'LineStyle','-','LineWidth',1);
    end
    
    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(-134,72);
end

