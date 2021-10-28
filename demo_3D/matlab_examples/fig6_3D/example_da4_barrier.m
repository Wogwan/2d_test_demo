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
c3 = (x1+0)^2+(x2-0)^2+(x3-6)^2-9;
c4 = (x1+0)^2+(x2+0)^2+(x3+5)^2-9;
    figure(1);clf;hold on;
    c1s = char(sdisplay(c1));
    c2s = char(sdisplay(c2));
    c3s = char(sdisplay(c3));
    c4s = char(sdisplay(c4));
    xran=[-8 8 -8 8 -8 8];
    smrplot(c1s,0,xran,[1000 300],'r-');
    smrplot(c2s,0,xran,[1000 300],'r-');
    smrplot(c3s,0,xran,[1000 300],'r-');
    smrplot(c4s,0,xran,[1000 300],'r-');
    refreshdata; drawnow;
% % largest level set
for ii = 1:10

% fix cc search for u, L1
[uu1 u1c u1v] = polynomial([x1 x2 x3], 1, 0);
[uu2 u2c u2v] = polynomial([x1 x2 x3], 1, 0);
[L1 L1c L1v] = polynomial([x1 x2 x3], 2, 0);

Vdot = jacobian(v, x1)*(x2-x3^2) + jacobian(v, x2)*(x3-x1^2+uu1) + jacobian(v, x3)*(-x1-2*x2-x3+x2^3+uu2);
F = [sos(L1), sos(-Vdot - L1*(cc - v)),max(abs([u1c; u2c]))<100];
[sol,vv,QQ] = solvesos(F,1,[],[L1c;u1c;u2c]);
%solL1 = clean(value(L1c)'*L1v,1e-6);
%solu = clean(value(uc)'*uv,1e-6);
solL1 = value(L1c)'*L1v;
solu1 = value(u1c)'*u1v;
solu2 = value(u2c)'*u2v;
sdisplay(solL1)
sdisplay(solu1)

% fix solu,solL1 search for vc
sdpvar vc
[uu1 u1c u1v] = polynomial([x1 x2 x3], 1, 0);
[uu2 u2c u2v] = polynomial([x1 x2 x3], 1, 0);
[L2 L2c L2v] = polynomial([x1 x2], 2, 0);
[L3 L3c L3v] = polynomial([x1 x2], 2, 0);
[L4 L4c L4v] = polynomial([x1 x2], 2, 0);
[L5 L5c L5v] = polynomial([x1 x2], 2, 0);

Vdot = jacobian(v, x1)*(x2-x3^2) + jacobian(v, x2)*(x3-x1^2+uu1) + jacobian(v, x3)*(-x1-2*x2-x3+x2^3+uu2);
F = [sos(L2),sos(L3), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c1*L2), sos(-(vc-v)+c2*L3), sos(-(vc-v)+c3*L4), sos(-(vc-v)+c4*L5),vc>0,max(abs([u1c; u2c]))<100];
[sol,vv,QQ] = solvesos(F,-vc,[],[u1c;u2c;L2c;L3c;L4c; L5c; vc]);
%solu = clean(value(uc)'*uv,1e-6);
solv = value(vc);

    cc = solv;
    % smr variables
    solh = value(vc) -v;
    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8];
    smrplot(hs,0,xran,[300 50],'y-');
    axis(xran)
    refreshdata; drawnow;
end

    figure(1);hold on;
    solh = value(vc) -v;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8];
    smrplot(hs,0,xran,[300 50],'b-');
    refreshdata; drawnow;

cc=13.0124;
% find barrier certificates
% init h with CLF
solh = cc - v;
    hinit = char(sdisplay(solh));
    xran=[-8 8 -8 8 -8 8];
    smrplot(hinit,0,xran,[300 50],'m-');
    axis(xran)
    refreshdata; drawnow;

for ii = 1:20
ii
sdpvar htol
% fix h search for u, L1, L2
[u1 u1c u1v] = polynomial([x1 x2 x3], 1, 1);
[u2 u2c u2v] = polynomial([x1 x2 x3], 1, 1);
[L1 L1c L1v] = polynomial([x1 x2 x3], 2, 0);
[L2 L2c L2v] = polynomial([x1 x2 x3], 2, 0);

hdot = jacobian(solh, x1)*(x2-x3^2) + jacobian(solh, x2)*(x3-x1^2+u1) + jacobian(solh, x3)*(-x1-2*x2-x3+x2^3+u2);
Vdot = jacobian(v, x1)*(x2-x3^2) + jacobian(v, x2)*(x3-x1^2+u1) + jacobian(v, x3)*(-x1-2*x2-x3+x2^3+u2);
F = [sos(L1), sos(L2), sos(-Vdot + L1*(-solh)), sos(hdot + gamma*solh + L2*(-solh)-htol), htol>=0, max(abs([u1c; u2c]))<=100];
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
[L5 L5c L5v] = polynomial([x1 x2 x3], 2, 0);
[L6 L6c L6v] = polynomial([x1 x2 x3], 2, 0);
hdot = jacobian(h, x1)*(x2-x3^2) + jacobian(h, x2)*(x3-x1^2+solu1) + jacobian(h, x3)*(-x1-2*x2-x3+x2^3+solu2);
Vdot = jacobian(v, x1)*(x2-x3^2) + jacobian(v, x2)*(x3-x1^2+solu1) + jacobian(v, x3)*(-x1-2*x2-x3+x2^3+solu2);
% F = [sos(L3), sos(L4),sos(-Vdot + solL1*(-h)), sos(hdot + gamma*h + solL2*(-h)), sos(-h+c1*L3), sos(-h+c2*L4), sos(-h+c3*L5), sos(-h+c4*L6)];
F = [sos(L3), sos(L4),sos(L5),sos(L6),sos(-Vdot + solL1*(-h)), sos(hdot + gamma*h + solL2*(-h)), sos(-h+c1*L3), sos(-h+c2*L4), sos(-h+c3*L5), sos(-h+c4*L6)];

[sol,vv,QQ] = solvesos(F,-hc(1)-hc(5)-hc(7)-hc(10),[],[L3c;L4c;L5c;L6c;hc]);
%solL3 = clean(value(L3c)'*L3v,1e-6);
%solL4 = clean(value(L4c)'*L4v,1e-6);
%solh = clean(value(hc)'*hv,1e-6);
temph = clean(value(hc)'*hv,1e-6);
if value(hc(1)) < 1
    break;
else
   solh = clean(value(hc)'*hv,1e-6);
end
sdisplay(solh)
    % smr variables
    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8 -8 8];
    smrplot(hs,0,xran,[1000 300],'m-');
    axis(xran)
    refreshdata; drawnow;
    view(-120,0)
end


    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8 -8 8];
    smrplot(hs,0,xran,[300 50],'b--');
    axis(xran)
    refreshdata; drawnow;

%% double check
hdot = jacobian(solh, x1)*(x2-x3^2) + jacobian(solh, x2)*(x3-x1^2+solu1) + jacobian(solh, x3)*(-x1-2*x2-x3+x2^3+solu2);

% solh:  114.3555451+1.4686*x1+7.2121*x2+19.8479*x3-24.5412*x3^2-14.7734*x1^2-26.0129*x1*x2-15.5440*x1*x3-28.3492*x2^2-27.5651*x2*x3
% solu1: -43.0365*x1-92.6311*x2-31.4178*x3
% solu2: -2.2003*x1-8.2927*x2-4.9452*x3
% optc: 13.0124
