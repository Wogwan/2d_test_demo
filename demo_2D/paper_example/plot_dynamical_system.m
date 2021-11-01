format long
hold on;
syms x1 x2;
dom = 4;
figure(15);clf;

hold on;
f = @(t,x)[
    -0.263205441705978*x(1)^2+1.861977417622602*x(1)*x(2)-1.956387277005004*x(2)^2-4207.835999213985*x(1)+5126.655463006833*x(2)-0.04136749982010112
    x(1)^2*x(2)+0.3898848123891965*x(1)^2-1.961131454038698*x(1)*x(2)+2.010549285003059*x(2)^2+2950.536581729438*x(1)-4081.561465971713*x(2)+0.03995820559230309
    ];
%% Plot ROA
% [ts1,ys1] = ode45(f,[0,20],[0.45;-0.51]);
% plot(ys1(:,1),ys1(:,2), 'r','linewidth',3), hold on;
%% Field [-2 2]
vectfield(f,-dom:0.4:dom,-dom:0.4:dom); hold on; % field [-5 5 -5 5]
for xd = -dom:dom/10:dom
    for yd = -dom:0.5:dom
        [ts,ys] = ode45(f,[0,20],[xd;yd]);
        h1 = plot(ys(:,1),ys(:,2), 'k');hold on;
    end
end
xlim([-dom dom]); ylim([-dom dom]);
plot(0,0,'c*','linewidth',2);hold on;
set(gca, 'LooseInset', [0,0,0,0]);
title('');