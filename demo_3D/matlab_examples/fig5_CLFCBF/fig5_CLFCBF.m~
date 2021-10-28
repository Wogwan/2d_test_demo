clc
clear
close all

% optimal c-v
v = 'x1^2+x1*x2+x2^2';
cinit = '0.1-x1^2-x1*x2-x2^2';
copt = '5.8628-x1^2-x1*x2-x2^2';
hopt = '0.518862456-0.0669*x1-0.1196*x2-0.0546*x1^2-0.0630*x1*x2-0.0294*x2^2';

% system equations
[xx yy] = meshgrid(-12:2:8);
xdot = yy; 
ydot = -xx + 3.245998424e-05 -0.1633*xx-2.0450*yy+0.0196*xx.^2+0.0685*xx.*yy+0.0689*yy.^2;
quiver(xx, yy, xdot, ydot,2,'m');
box off;

xran=[-6 8 -12 8];
hold on;
smrplot(copt,0,xran,[400 100],'g-');
hs = findobj('Type', 'line');
%set(hs, 'Linewidth', 2);
smrplot(hopt,0,xran,[300 50],'b--');
hs = findobj('Type', 'line');
set(hs, 'Linewidth', 3);
legend(hs([1,end]), {'$h^*(x)=0$','$V(x)=c^*$'}, 'Interpreter','latex');
axis(xran);


% safety constraints
sdpvar x1 x2
c1 = (x1-3)^2+(x2-1)^2-1;
c2 = (x1+3)^2+(x2+4)^2-1;
c3 = (x1+4)^2+(x2-5)^2-1;
c4 = (x1-2)^2+(x2+6)^2-1;
hold on;
c1s = char(sdisplay(c1));
c2s = char(sdisplay(c2));
c3s = char(sdisplay(c3));
c4s = char(sdisplay(c4));
smrplot(c1s,0,xran,[300 50],'r-');
smrplot(c2s,0,xran,[300 50],'r-');
smrplot(c3s,0,xran,[300 50],'r-');
%smrplot(c4s,0,xran,[300 50],'r-');
%smrplot(cinit,0,xran,[300 50],'c-');





xlabel('$x_1$', 'Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
axis(xran)
set(gcf, 'color', [1 1 1]);
set(gca, 'LooseInset', [0,0,0,0]);
x0=1;y0=0;width=3.5;height=3.5*0.75;
set(gcf,'units','inches','position',[x0,y0,width,height]) % set window size

