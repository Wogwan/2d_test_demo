clear; 
format long
hold on;
syms x1 x2;
dom = 4;
figure(16);clf;

hold on;
% f = @(t,x)[
%     -x(1)+x(2)
%     x(1)^2*x(2)+1-sqrt(sqrt((exp(x(1))*cos(x(1)))^2))
%     ];

f = @(t,x)[
    x(2)-x(1)
    -0.48494305021249242693137659898639*x(1)^4+0.4476352069538196976061783516343*x(1)^3+x(1)^2*x(2)+0.35813066874120685900706462234666*x(1)^2-0.47096752768685094803213786462948*x(1)+0.03729676680157056195552556232542*x(2)^4+0.01893812978180069162004173222158*x(2)^3+0.13680422363948308017711497086566*x(2)^2+0.043563863537901807709840085180986*x(2)-0.33514268478706185638849035512976
    ];

%% Plot ROA
[ts1,ys1] = ode45(f,[0,20],[0.1;-0.1]);
plot(ys1(:,1),ys1(:,2), 'r','linewidth',3), hold on;
%% Field [-2 2]
vectfield(f,-dom:0.4:dom,-dom:0.4:dom); hold on; % field [-5 5 -5 5]
for xd = -dom:dom/20:dom
    for yd = -dom:dom/20:dom
        [ts,ys] = ode45(f,[0,20],[xd;yd]);
        h1 = plot(ys(:,1),ys(:,2), 'k');hold on;
    end
end
xlim([-dom dom]); ylim([-dom dom]);
plot(0,0,'c*','linewidth',2);hold on;
set(gca, 'LooseInset', [0,0,0,0]);
title('');