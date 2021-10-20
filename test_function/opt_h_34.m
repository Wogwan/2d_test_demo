clc
clear
close all

sdpvar x1 x2 u
f = [-x2-3/2*x1^2-1/2*x1^3; x1 - u];
v = x1^2+x2^2+1*x1*x2;
e = [0;0];
c = 3.00000001711;
gamma = 1;

% init h with CLF
solh = c - v;
% safety constraints
c1 = (x1-3)^2+(x2-1)^2-1;
c2 = (x1+3)^2+(x2+4)^2-1;
c3 = (x1+4)^2+(x2-5)^2-1;
c4 = (x1-2)^2+(x2+6)^2-1;
c5 = -x2+2;
c6 =(x1+0.5)^2+(x2+3.25)^2-1;
    figure(2);clf;hold on;
    c1s = char(sdisplay(c1));
    c2s = char(sdisplay(c2));
    c3s = char(sdisplay(c3));
    c4s = char(sdisplay(c4));
    c5s = char(sdisplay(c5));
    c6s = char(sdisplay(c6));
    hinit = char(sdisplay(solh));
    dom = 20;
    xran=[-dom dom -dom dom];
%     smrplot(c1s,0,xran,[300 50],'r-');
%     smrplot(c2s,0,xran,[300 50],'r-');
%     smrplot(c3s,0,xran,[300 50],'r-');
    %smrplot(c4s,0,xran,[300 50],'r-');
    smrplot(c5s,0,xran,[1000 400],'r-');
    smrplot(c6s,0,xran,[1000 400],'r-');
    smrplot(hinit,0,xran,[1000 400],'c-');
    axis(xran)


while 1
sdpvar htol
% fix h search for u, L1, L2
%[u, L1, L2] = sossearch_u(solh)
[u uc uv] = polynomial([x1 x2], 2, 0);
[L1 L1c L1v] = polynomial([x1 x2], 2, 0);
[L2 L2c L2v] = polynomial([x1 x2], 2, 0);

hdot = jacobian(solh, x1)*x2 + jacobian(solh, x2)*(-x1 + u);
Vdot = jacobian(v, x1)*x2 + jacobian(v, x2)*(-x1 + u);
F = [sos(L1), sos(L2), sos(-Vdot + L1*(-solh)), sos(hdot + gamma*solh + L2*(-solh)-htol), htol>=0];
[sol,vv,QQ] = solvesos(F,-htol,[],[L1c;L2c;uc]);
%solL1 = clean(value(L1c)'*L1v,1e-6);
%solL2 = clean(value(L2c)'*L2v,1e-6);
%solu = clean(value(uc)'*uv,1e-6);
solL1 = value(L1c)'*L1v;
solL2 = value(L2c)'*L2v;
solu = value(uc)'*uv;
sdisplay(solL1)
sdisplay(solL2)
sdisplay(solu)
%3.245998424e-05-0.1633*x1-2.0450*x2+0.0196*x1^2+0.0685*x1*x2+0.0689*x2^2

% fix u, L1, L2, search for h, L3, L4
%[h] = sossearch_h(solu, solL1, solL2)
[h hc hv] = polynomial([x1 x2], 2, 0);
% [L3 L3c L3v] = polynomial([x1 x2], 2, 0);
% [L4 L4c L4v] = polynomial([x1 x2], 2, 0);
[L5 L5c L5v] = polynomial([x1 x2], 2, 0);
[L6 L6c L6v] = polynomial([x1 x2], 2, 0);

hdot = jacobian(h, x1)*x2 + jacobian(h, x2)*(-x1 + solu);
Vdot = jacobian(v, x1)*x2 + jacobian(v, x2)*(-x1 + solu);
% F = [sos(L3), sos(L4), sos(L5), sos(-Vdot + solL1*(-h)), sos(hdot + gamma*h + solL2*(-h)), sos(-h+c1*L3), sos(-h+c2*L4), sos(-h+c3*L5)];
% [sol,vv,QQ] = solvesos(F,-hc(1)-hc(4)-hc(6),[],[L3c;L4c;L5c; hc]);
F = [sos(L5), sos(L6), sos(-Vdot + solL1*(-h)), sos(hdot + gamma*h + solL2*(-h)), sos(-h+c5*L5), sos(-h+c6*L6)];
[sol,vv,QQ] = solvesos(F,-hc(1)-hc(4)-hc(6),[],[L5c;L6c;hc]);
%-hc(1)-hc(4)-hc(6)-hc(11)-hc(13)-hc(15)
%solL3 = clean(value(L3c)'*L3v,1e-6);
%solL4 = clean(value(L4c)'*L4v,1e-6);
%solh = clean(value(hc)'*hv,1e-6);
% solL3 = value(L3c)'*L3v;
% solL4 = value(L4c)'*L4v;
solL5 = value(L5c)'*L5v;
solL6 = value(L6c)'*L6v;
solh = value(hc)'*hv;
% sdisplay(solL3)
% sdisplay(solL4)
sdisplay(solh)
% 0.518862456-0.0669*x1-0.1196*x2-0.0546*x1^2-0.0630*x1*x2-0.0294*x2^2
    % smr variables
    figure(2);hold on;
    hs = char(sdisplay(solh));
    xran=[-dom dom -dom dom];
    smrplot(hs,0,xran,[1000 400],'g-');
    axis(xran)
% solu: 3.245998424e-05-0.1633*x1-2.0450*x2+0.0196*x1^2+0.0685*x1*x2+0.0689*x2^2
end


    figure(2);hold on;
%     hs = char(sdisplay(solh));
%     xran=[-dom dom -dom dom];
    smrplot(hs,0,xran,[1000 400],'b-');
    axis(xran)
