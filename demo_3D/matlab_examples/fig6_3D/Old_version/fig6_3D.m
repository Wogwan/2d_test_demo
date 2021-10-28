clc
close all
clear all

sdpvar x1 x2 x3 u1 u2
f = [x2-x3^2; x3-x1^2+u1; -x1-2*x2-x3+x2^3+u2];
v = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
copt = 28.8735;
vc = copt - v;
gamma = 2;
hopt = 36.66082419+0.7231*x1+0.9507*x2+0.6429*x3-5.4803*x3^2 -3.5685*x1^2-10.8668*x1*x2-3.9342*x1*x3-11.1044*x2^2-10.0481*x2*x3;

vcs = char(sdisplay(vc));
xran=[-7 7 -5 5 -4 4];
smrplot(vcs,0,xran,[100 50 40],'k-');
hold on;

hs = char(sdisplay(hopt));
smrplot(hs,0,xran,[100 30 40],'b-');

% safety constraints
c1 = (x1-2)^2+(x2-1)^2+(x3-2)^2-1;
c2 = (x1+1)^2+(x2+2)^2+(x3+1)^2-1;
hold on;
c1s = char(sdisplay(c1));
c2s = char(sdisplay(c2));
smrplot(c1s,0,xran,[100 50 40],'r-');
smrplot(c2s,0,xran,[100 50 40],'r-');
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');
%set(gcf, 'color', [1 1 1]);
%set(gca, 'LooseInset', [0,0,0,0]);
x0=1;y0=0;width=5;height=5*0.75;
set(gcf,'units','inches','position',[x0,y0,width,height]) % set window size
