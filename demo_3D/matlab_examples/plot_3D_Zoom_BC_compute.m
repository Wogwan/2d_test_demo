% clc;
% load BC_test_Lya;
% load BC_test_Lya2;
% load test_1125_1_3D;
% load test_1125_3_3D;
% load figure_3D;
clear;
load test_1125_night_3D;
Barrier = Compute{7};
pvar x1 x2 x3;
dom = 10; domain = [-dom dom -dom dom -dom dom];
C1 = (x1-3)^2+(x2-0)^2+(x3-0)^2-2;
C2 = (x1-0)^2+(x2+3)^2+(x3+0)^2-2;
C3 = (x1-3)^2+(x2-2)^2+(x3-2)^2-2;
C4 = (x1+3)^2+(x2-3)^2+(x3+2)^2-2;
V = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
C0 = 6.323809779487249;
k = ['r','g','m','c','k','y'];
%%
% for iter = 5:length(A)
% for iter = 30
%     B = Barrier(iter,1);
%     B = A(iter,4);
iter = 1;
Ba = Barrier{iter};
for iter = 1:length(Ba)
    B = Ba(iter,1);
    %     B = 5.107896386952879e-05*x1^4-1.225600219656938e-06*x1^3*x2+1.276555528738019e-05*x1^3*x3+0.0001329541596461657*x1^2*x2^2-2.877737078915467e-06*x1^2*x2*x3+5.636341555815255e-05*x1^2*x3^2-1.334560486464528e-06*x1*x2^3-1.515539026477243e-05*x1*x2^2*x3-1.259767196325635e-08*x1*x2*x3^2+8.043987202842057e-06*x1*x3^3-1.856106327054302e-05*x2^4-1.394556287044935e-06*x2^3*x3+6.995213986226199e-05*x2^2*x3^2+1.569022379844585e-06*x2*x3^3-8.714438896588289e-06*x3^4-3.690472179505634e-05*x1^3+4.615166513438906e-06*x1^2*x2-1.166030133789836e-06*x1^2*x3+0.0001141265996146262*x1*x2^2+4.245453997185847e-05*x1*x2*x3+2.718602832031601e-05*x1*x3^2+2.783412212246791e-06*x2^3+4.704749994552823e-05*x2^2*x3+1.543839088154108e-05*x2*x3^2-1.495660096986019e-05*x3^3-0.003041091903546785*x1^2+7.130821578404228e-05*x1*x2-0.0001335869126843515*x1*x3-0.001659588818064906*x2^2-0.0001774918790187125*x2*x3-0.0001076361419015228*x3^2-0.0005011300570146603*x1-0.0001011785249374663*x2+5.49385344531481e-06*x3+0.01000507213238431;
%     figure(99999);
    figure(iter);clf;hold on;
    phB= patch(pcontour3(B,0,domain,'B'));
    phV0= patch(pcontour3(V,double(C0),domain,'c')); set(phV0,'EdgeAlpha',1, 'FaceColor', 'none', 'EdgeColor', 'b' );
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

