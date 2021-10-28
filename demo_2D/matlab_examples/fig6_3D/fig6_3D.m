clc
close all
clear all

sdpvar x1 x2 x3 u1 u2
f = [x2-x3^2; x3-x1^2+u1; -x1-2*x2-x3+x2^3+u2];
v = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
copt = 13.0124;
vc = copt - v;
gamma = 2;
hopt = 114.3555451+1.4686*x1+7.2121*x2+19.8479*x3-24.5412*x3^2-14.7734*x1^2-26.0129*x1*x2-15.5440*x1*x3-28.3492*x2^2-27.5651*x2*x3;

vcs = char(sdisplay(vc));
xran=[-7 7 -5 5 -4 4];
smrplot(vcs,0,xran,[100 50 40],'k-');
hold on;
% solu1: -43.0365*x1-92.6311*x2-31.4178*x3
% solu2: -2.2003*x1-8.2927*x2-4.9452*x3


hs = char(sdisplay(hopt));
smrplot(hs,0,xran,[100 30 40],'b-');

% safety constraints
c1 = (x1-2)^2+(x2-1)^2+(x3-2)^2-1;
c2 = (x1+1)^2+(x2+2)^2+(x3+1)^2-1;
c3 = (x1+0)^2+(x2-0)^2+(x3-6)^2-9;
c4 = (x1+0)^2+(x2+0)^2+(x3+5)^2-9;
hold on;
c1s = char(sdisplay(c1));
c2s = char(sdisplay(c2));
c3s = char(sdisplay(c3));
c4s = char(sdisplay(c4));
smrplot(c1s,0,xran,[100 50 40],'r-');
smrplot(c2s,0,xran,[100 50 40],'r-');
smrplot(c3s,0,xran,[100 50 40],'r-');
smrplot(c4s,0,xran,[100 50 40],'r-');
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');
set(gcf, 'color', [1 1 1]);
%set(gca, 'LooseInset', [0,0,0,0]);
x0=1;y0=0;width=5;height=5*0.75;
set(gcf,'units','inches','position',[x0,y0,width,height]) % set window size
