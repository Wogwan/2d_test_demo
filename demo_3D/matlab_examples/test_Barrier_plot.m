% figure_id = 22;
% figure(figure_id);clf;hold on;
pvar x1 x2 x3
dom = 10;
domain = [-dom dom -dom dom -dom dom];
for i = 1:length(Barrier)
    % for i = 1:length(A)
    i
    figure(i+300);clf;hold on;
    V = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
    C0 = 9.681159722441127;
    us4 = patch(pcontour3(V,C0,domain,'r'));
    set(us4, 'EdgeAlpha',0.4,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','--','LineWidth',1 ); hold on;
    inH = patch(pcontour3(Barrier(i),0,domain,'k'));
    %     inH = patch(pcontour3(A(i,4),0,domain,'k'));
    set(inH, 'EdgeAlpha',0.6,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','-','LineWidth',0.8 ); hold on;
    
    
    view(-185,5);hold on;    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]);
end


% % Original
% inH = patch(pcontour3(A(i+9,4),0,domain,'k'));
% set(inH, 'EdgeAlpha',0.6,'FaceColor', 'none', 'EdgeColor', 'k','LineStyle','-','LineWidth',0.8 ); hold on;
% %% Another
% us4 = patch(pcontour3(Barrier(i),0,domain,'k'));
% set(us4, 'EdgeAlpha',0.4,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','--','LineWidth',1 ); hold on;
