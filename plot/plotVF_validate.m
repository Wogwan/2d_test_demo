function plotVF_validate

% hold on
format long
hold on;
syms x1 x2;
sample_time = 40;
dom = 10;

%%
% GP system
figure(11);hold on;
f = @(t,x) [
-x(2)-3/2*x(1)^2-1/2*x(1)^3
x(1)];

%%
% % GP system
% figure(22);clf;
% f = @(t,x) [
% x(2)-x(1)
% x(1)^2*x(2)-0.000007031474139-0.500009736044418*x(1)+0.126057255453549*x(1)^2+0.103008552092872*x(1)^3-0.168853172050168 *x(1)^4+0.264932312148793*x(1)^5;
% ];

%% Plot ROA
% [ts1,ys1] = ode45(f,[0,20],[0.45;-0.51]);
% plot(ys1(:,1),ys1(:,2), 'r','linewidth',3), hold on;


%% Field [-2 2]
% vectfield(f,-2:0.4:2,-2:0.4:2); hold on; % field [-5 5 -5 5]
vectfield(f,-dom:dom/10:dom,-dom:dom/10:dom); hold on; % field [-5 5 -5 5]


for xd = -dom:dom/10:dom
    for yd = -dom:dom/10:dom
        [ts,ys] = ode45(f,[0,sample_time],[xd;yd]);
        h1 = plot(ys(:,1),ys(:,2), 'k');hold on;
    end
end

xlim([-dom dom]); ylim([-dom dom]);
% plot(0,0,'co','linewidth',3);hold on;

set(gca, 'LooseInset', [0,0,0,0]);
title('');

end