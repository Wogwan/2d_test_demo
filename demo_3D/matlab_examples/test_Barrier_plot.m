figure_id = 22;
figure(figure_id);clf;hold on;
dom = 15;
domain = [-dom dom -dom dom -dom dom];
for i = 1:length(Barrier)
    i
    figure(i);clf;hold on;
    inH = patch(pcontour3(Barrier(i),0,domain,'k'));
    set(inH, 'EdgeAlpha',1,'FaceColor', 'none', 'EdgeColor', 'm','LineStyle','-','LineWidth',2 ); hold on;
    view(-30,20);hold on;    xlim([-dom dom]); ylim([-dom dom]); zlim([-dom dom]); view(-30,20);
end

% 8?