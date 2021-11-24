clear;close all;clc;
% load BC_test_Lya;
load BC_test_Lya2;
dom = 10; domain = [-dom dom -dom dom -dom dom];

k = ['r','g','b','m','c','k','y'];
for iter = 1:length(A)
    B = A(iter,5);
    figure(iter);clf;hold on;
    phB= patch(pcontour3(B,0,domain,'B'));
    
    if mod(iter,7) == 0
        set(phB,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', k(7),'LineStyle','-','LineWidth',1);
    else
        set(phB,'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', k(mod(iter,7)),'LineStyle','-','LineWidth',1);
    end
    
    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(-30,20);
    pause
end

