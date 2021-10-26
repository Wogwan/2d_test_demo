clc
clear
close all

sdpvar x1 x2 u
f = [x2; -x1 + u];
v = x1^2+x1*x2+x2^2;
e = [0;0];
cc = 0.1;
gamma = 1;

% init h with CLF
solh = cc - v;
% safety constraints
c1 = (x1-3)^2+(x2-1)^2-1;
c2 = (x1+3)^2+(x2+4)^2-1;
c3 = (x1+4)^2+(x2-5)^2-1;
c4 = (x1-2)^2+(x2+6)^2-1;
    figure(1);clf;hold on;
    c1s = char(sdisplay(c1));
    c2s = char(sdisplay(c2));
    c3s = char(sdisplay(c3));
    c4s = char(sdisplay(c4));
    hinit = char(sdisplay(solh));
    xran=[-8 8 -8 8];
    smrplot(c1s,0,xran,[300 50],'r-');
    smrplot(c2s,0,xran,[300 50],'b-');
    smrplot(c3s,0,xran,[300 50],'g-');
    %smrplot(c4s,0,xran,[300 50],'r-');
    smrplot(hinit,0,xran,[300 50],'c-');
    axis(xran)


for ii = 1:10

% fix cc search for u, L1
[u uc uv] = polynomial([x1 x2], 2, 0);
[L1 L1c L1v] = polynomial([x1 x2], 2, 0);

Vdot = jacobian(v, x1)*x2 + jacobian(v, x2)*(-x1 + u);
F = [sos(L1), sos(-Vdot - L1*(cc - v))];
[sol,vv,QQ] = solvesos(F,1,[],[L1c;uc]);
%solL1 = clean(value(L1c)'*L1v,1e-6);
%solu = clean(value(uc)'*uv,1e-6);
solL1 = value(L1c)'*L1v;
solu = value(uc)'*uv;
sdisplay(solL1)
sdisplay(solu)

% fix solu,solL1 search for vc
sdpvar vc
[u uc uv] = polynomial([x1 x2], 2, 0);
[L2 L2c L2v] = polynomial([x1 x2], 2, 0);
[L3 L3c L3v] = polynomial([x1 x2], 2, 0);
[L4 L4c L4v] = polynomial([x1 x2], 2, 0);

Vdot = jacobian(v, x1)*x2 + jacobian(v, x2)*(-x1 + u);
%F = [sos(L2),sos(L3),sos(L4), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c1*L2), sos(-(vc-v)+c2*L3), sos(-(vc-v)+c2*L4), vc>0];
%%
% F = [sos(L2),sos(L3),sos(L4), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c1*L2), sos(-(vc-v)+c2*L3), sos(-(vc-v)+c2*L4), vc>0];
% [sol,vv,QQ] = solvesos(F,-vc,[],[uc;L2c; L3c; L4c;vc]);
%%
F = [sos(L2), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c3*L2), vc>0];
[sol,vv,QQ] = solvesos(F,-vc,[],[uc;L2c;vc]);

%solu = clean(value(uc)'*uv,1e-6);
solv = value(vc);

    cc = solv;
    % smr variables
    solh = value(vc) -v;
    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8];
    smrplot(hs,0,xran,[300 50],'g-');
    axis(xran)

end


    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8];
    figure(1);hold on;
    smrplot(hs,0,xran,[300 50],'b-');
    axis(xran)
