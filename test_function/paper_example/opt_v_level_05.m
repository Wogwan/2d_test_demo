clear
close all

sdpvar x1 x2 u
%%
% f = [x2; -x1 + u];
% v = x1^2+x1*x2+x2^2;
%%
% f = [-x2-3/2*x1^2-1/2*x1^3; x1-u];
% % v = 1*x1^2+4*x2^2+1*x1*x2;
% v = 4*x1^2+8*x2^2-4*x1*x2+2*x1+2*x2;
%%
f = [0.1*x1^2+1*x2; 0.1*x1*x2-0.2*x1+(1+x1^2)*u];
v = 1*x1^2+1*x2^2+1*x1*x2;

%%
e = [0;0];
cc = 2;
gamma = 0;
dom = 20;

% init h with CLF
solh = cc - v;
% safety constraints
c1 = (x1-3)^2+(x2-1)^2-1;
c2 = (x1+3)^2+(x2+4)^2-1;
c3 = (x1+4)^2+(x2-5)^2-1;
c4 = (x1-2)^2+(x2-6)^2-1;
c5 = -x2+2;
c6 =(x1+0.5)^2+(x2+3.25)^2-1;
    figure(1);clf;hold on;
%     c1s = char(sdisplay(c 1));
%     c2s = char(sdisplay(c2));
%     c3s = char(sdisplay(c3));
    c4s = char(sdisplay(c4));
%     c5s = char(sdisplay(c5));
    c6s = char(sdisplay(c6));
    hinit = char(sdisplay(solh));
    xran=[-dom dom -dom dom];
%     smrplot(c1s,0,xran,[300 50],'r-');
%     smrplot(c2s,0,xran,[300 50],'r-');
%     smrplot(c3s,0,xran,[300 50],'r-');
    smrplot(c4s,0,xran,[300 50],'r-');
%     smrplot(c5s,0,xran,[1000 300],'r-');
    smrplot(c6s,0,xran,[1000 300],'r-');
    smrplot(hinit,0,xran,[1000 300],'c-');
    axis(xran)


% for ii = 1:10
while 1

% fix cc search for u, L1
[u uc uv] = polynomial([x1 x2], 2, 0);
[L1 L1c L1v] = polynomial([x1 x2], 2, 0);

% Vdot = jacobian(v, x1)*f(1) + jacobian(v, x2)*f(2);
Vdot = jacobian(v, x1)*(0.1*x1^2+1*x2)+jacobian(v,x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*u);
F = [sos(L1), sos(-Vdot - L1*(cc - v))];
[sol,vv,QQ] = solvesos(F,1,[],[L1c;uc]);
% solL1 = clean(value(L1c)'*L1v,1e-6);
% solu = clean(value(uc)'*uv,1e-6);
solL1 = value(L1c)'*L1v;
solu = value(uc)'*uv;
sdisplay(solL1)
sdisplay(solu)

% fix solu,solL1 search for vc
sdpvar vc
% [u uc uv] = polynomial([x1 x2], 2, 0);
% [L2 L2c L2v] = polynomial([x1 x2], 2, 0);
% [L3 L3c L3v] = polynomial([x1 x2], 2, 0);
[L4 L4c L4v] = polynomial([x1 x2], 2, 0);
% [L5 L5c L5v] = polynomial([x1 x2], 2, 0);
[L6 L6c L6v] = polynomial([x1 x2], 2, 0);

Vdot = jacobian(v, x1)*(0.1*x1^2+1*x2)+jacobian(v,x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*solu);
% Vdot = jacobian(v, x1)*f(1) + jacobian(v, x2)*f(2);
% Vdot = jacobian(v, x1)*x2 + jacobian(v, x2)*(-x1 + solu);
%F = [sos(L2),sos(L3),sos(L4), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c1*L2), sos(-(vc-v)+c2*L3), sos(-(vc-v)+c2*L4), vc>0];
% F = [sos(L2),sos(L3),sos(L4), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c1*L2), sos(-(vc-v)+c2*L3), sos(-(vc-v)+c2*L4), vc>0];
F = [sos(L4),sos(L6), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c4*L4), sos(-(vc-v)+c6*L6), vc>0];
[sol,vv,QQ] = solvesos(F,-vc,[],[uc;L4c; L6c; vc]);
%solu = clean(value(uc)'*uv,1e-6);
solv = value(vc);

    cc = solv;
    % smr variables
    solh = value(vc) -v;
    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-dom dom -dom dom];
    smrplot(hs,0,xran,[1000 300],'g-');
    axis(xran)

end


    figure(1);hold on;
    hs = char(sdisplay(solh))
    xran=[-dom dom -dom dom];
    smrplot(hs,0,xran,[1000 300],'b-');
    axis(xran)
