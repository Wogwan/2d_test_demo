clc
close all
clear all

%Init system
sdpvar x1 x2 x3 u1 u2
f = [x2-x3^2; x3-x1^2+u1; -x1-2*x2-x3+x2^3+u2];
v = 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2;
cc = 1;
gamma = 2;

% safety constraints
c1 = (x1-2)^2+(x2-1)^2+(x3-2)^2-1;
c2 = (x1+1)^2+(x2+2)^2+(x3+1)^2-1;
    hold on;
    c1s = char(sdisplay(c1));
    c2s = char(sdisplay(c2));
    xran=[-8 8 -8 8 -8 8];
    smrplot(c1s,0,xran,[300 50],'r-');
    smrplot(c2s,0,xran,[300 50],'r-');



% % largest level set
% for ii = 1:20

% % fix cc search for u, L1
% [u1 u1c u1v] = polynomial([x1 x2 x3], 1, 0);
% [u2 u2c u2v] = polynomial([x1 x2 x3], 1, 0);
% [L1 L1c L1v] = polynomial([x1 x2 x3], 2, 0);

% Vdot = jacobian(v, x1)*(x2-x3^2) + jacobian(v, x2)*(x3-x1^2+u1) + jacobian(v, x3)*(-x1-2*x2-x3+x2^3+u2);
% F = [sos(L1), sos(-Vdot - L1*(cc - v))];
% [sol,vv,QQ] = solvesos(F,1,[],[L1c;u1c;u2c]);
% %solL1 = clean(value(L1c)'*L1v,1e-6);
% %solu = clean(value(uc)'*uv,1e-6);
% solL1 = value(L1c)'*L1v;
% solu1 = value(u1c)'*u1v;
% solu2 = value(u2c)'*u2v;
% sdisplay(solL1)
% sdisplay(solu1)

% % fix solu,solL1 search for vc
% sdpvar vc
% [u1 u1c u1v] = polynomial([x1 x2 x3], 1, 0);
% [u2 u2c u2v] = polynomial([x1 x2 x3], 1, 0);
% [L2 L2c L2v] = polynomial([x1 x2], 2, 0);
% [L3 L3c L3v] = polynomial([x1 x2], 2, 0);

% Vdot = jacobian(v, x1)*(x2-x3^2) + jacobian(v, x2)*(x3-x1^2+u1) + jacobian(v, x3)*(-x1-2*x2-x3+x2^3+u2);
% F = [sos(L2),sos(L3), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c1*L2), sos(-(vc-v)+c2*L3),vc>0];
% [sol,vv,QQ] = solvesos(F,-vc,[],[u1c;u2c;L2c;L3c;vc]);
% %solu = clean(value(uc)'*uv,1e-6);
% solv = value(vc);

%     cc = solv;
%     % smr variables
%     solh = value(vc) -v;
%     hold on;
%     hs = char(sdisplay(solh));
%     xran=[-8 8 -8 8];
%     smrplot(hs,0,xran,[300 50],'r-');
%     axis(xran)

% end

    hold on;
    solh = value(28.8735) -v;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8];
    smrplot(hs,0,xran,[300 50],'r-');


cc=20;
% find barrier certificates
% init h with CLF
solh = cc - v;
    hinit = char(sdisplay(solh));
    xran=[-8 8 -8 8 -8 8];
    smrplot(hinit,0,xran,[300 50],'c-');
    axis(xran)


for ii = 1:3
ii
sdpvar htol
% fix h search for u, L1, L2
[u1 u1c u1v] = polynomial([x1 x2 x3], 1, 0);
[u2 u2c u2v] = polynomial([x1 x2 x3], 1, 0);
[L1 L1c L1v] = polynomial([x1 x2 x3], 2, 0);
[L2 L2c L2v] = polynomial([x1 x2 x3], 2, 0);

hdot = jacobian(solh, x1)*(x2-x3^2) + jacobian(solh, x2)*(x3-x1^2+u1) + jacobian(solh, x3)*(-x1-2*x2-x3+x2^3+u2);
Vdot = jacobian(v, x1)*(x2-x3^2) + jacobian(v, x2)*(x3-x1^2+u1) + jacobian(v, x3)*(-x1-2*x2-x3+x2^3+u2);
F = [sos(L1), sos(L2), sos(-Vdot + L1*(-solh)), sos(hdot + gamma*solh + L2*(-solh)-htol), htol>=0];
[sol,vv,QQ] = solvesos(F,-htol,[],[L1c;L2c;u1c;u2c]);
%solL1 = clean(value(L1c)'*L1v,1e-6);
%solL2 = clean(value(L2c)'*L2v,1e-6);
%solu = clean(value(uc)'*uv,1e-6);
solL1 = value(L1c)'*L1v;
solL2 = value(L2c)'*L2v;
solu1 = value(u1c)'*u1v;
solu2 = value(u2c)'*u2v;
%sdisplay(solL1)
%sdisplay(solL2)
%sdisplay(solu1)
%sdisplay(solu2)



% fix u, L1, L2, search for h, L3, L4
[h hc hv] = polynomial([x1 x2 x3], 2, 0);
[L3 L3c L3v] = polynomial([x1 x2 x3], 2, 0);
[L4 L4c L4v] = polynomial([x1 x2 x3], 2, 0);
hdot = jacobian(h, x1)*(x2-x3^2) + jacobian(h, x2)*(x3-x1^2+solu1) + jacobian(h, x3)*(-x1-2*x2-x3+x2^3+solu2);
Vdot = jacobian(v, x1)*(x2-x3^2) + jacobian(v, x2)*(x3-x1^2+solu1) + jacobian(v, x3)*(-x1-2*x2-x3+x2^3+solu2);
F = [sos(L3), sos(L4),sos(-Vdot + solL1*(-h)), sos(hdot + gamma*h + solL2*(-h)), sos(-h+c1*L3), sos(-h+c2*L4)];
[sol,vv,QQ] = solvesos(F,-hc(1)-hc(5)-hc(7)-hc(10),[],[L3c;L4c;hc]);
%solL3 = clean(value(L3c)'*L3v,1e-6);
%solL4 = clean(value(L4c)'*L4v,1e-6);
%solh = clean(value(hc)'*hv,1e-6);
temph = clean(value(hc)'*hv,1e-6);
if value(temph) == 0
    break;
else
   solh = clean(value(hc)'*hv,1e-6);
end
sdisplay(solh)
    % smr variables
    hold on;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8 -8 8];
    smrplot(hs,0,xran,[300 50],'g--');
    axis(xran)

end

%optimal h
%36.66082419+0.7231*x1+0.9507*x2+0.6429*x3-5.4803*x3^2-3.5685*x1^2-10.8668*x1*x2-3.9342*x1*x3-11.1044*x2^2-10.0481*x2*x3
% u1: -0.607104273-2.4861e+03*x1-3.3909e+03*x2+3.8618e+03*x3
% u2: 2.133083556+2.9243e+03*x1+3.4286e+03*x2-6.5492e+03*x3

    hold on;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8 -8 8];
    smrplot(hs,0,xran,[300 50],'b--');
    axis(xran)


