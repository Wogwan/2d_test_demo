clear; close all;
format long
hold on;
syms x1 x2;
dom = 4;
figure(15);clf;

hold on;
f = @(t,x)[
    -x(1)+x(2)
    x(1)^2*x(2)+1-sqrt(sqrt((exp(x(1))*cos(x(1)))^2))
    ];

% f = @(t,x)[
% x(2)-x(1)
% -0.16360815023548774815864703668922*x(1)^4-0.36060788497674955976890487363562*x(1)^3+x(1)^2*x(2)-0.058756810537387405002363038875046*x(1)^2-0.55810620811608024904870717364247*x(1)+0.78211290467871807940980488638161*x(2)^4+0.27124198816108746612485447258223*x(2)^3-0.015358122714423431964814170669342*x(2)^2+0.026176523031019402476538004975737*x(2)-0.33752390229169604296544093813282
%     ];

%% Plot ROA
[ts1,ys1] = ode45(f,[0,20],[0.1;-0.1]);
plot(ys1(:,1),ys1(:,2), 'r','linewidth',3), hold on;
%% Field [-2 2]
vectfield(f,-dom:0.4:dom,-dom:0.4:dom); hold on; % field [-5 5 -5 5]
for xd = -dom:dom/12:dom
    for yd = -dom:0.5:dom
        [ts,ys] = ode45(f,[0,20],[xd;yd]);
        h1 = plot(ys(:,1),ys(:,2), 'k');hold on;
    end
end
xlim([-dom dom]); ylim([-dom dom]);
plot(0,0,'c*','linewidth',2);hold on;
set(gca, 'LooseInset', [0,0,0,0]);
title('');