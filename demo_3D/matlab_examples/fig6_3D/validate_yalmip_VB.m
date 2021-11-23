clear all
close all
%Init system
sdpvar x1 x2 x3 u1 u2
f = [-0.16211179709695037165964324641892*x1^4+0.49031898059485488031675707268867*x1^3-0.80741059860165544525334054266171*x1^2-0.67918469453797989758096302163419*x1-0.000031796033448678815965058458425929*x2^4+0.00049811603934539429219124917480599*x2^3-0.017884574789259536503616132563366*x2^2+0.059003602596920091960530641017613*x2+0.14736298331076760903535216584714*x3^4-0.22234685217253638556123007674614*x3^3+0.14008297734444075111071015271591*x3^2+0.34076230832989312657943514750514*x3-1.4865840251194628709408125437986
    -x2*x1^3-x2
    1.876321954439673866943394386908*x1^4-1.8803947547684580765547934788628*x1^3-x1^2*x3-0.3521313244010295662178577913437*x1^2-0.58570314276461321600919518459705*x1+0.0008606977977094753011130801034767*x2^4-0.0047246616639717428295930368165045*x2^3+0.0073970075005580443808228530144788*x2^2-0.39277632595769917944750204696902*x2-0.63399852647351906398398568853736*x3^4-2.3526877462367330462456038731034*x3^3-0.18912599086086179234200699283974*x3^2-1.2624646738177363047839207865763*x3+1.2361363889678920191528277428006
    ];
v = 1*x1^4+1*x2^4+1*x3^4+1*x1^2*x2^2+1*x3^2*x2^2+1*x1^2*x3^2;
cc = 0.1;
gamma = 0;
dom = 10;

% safety constraints
% c1 = (x1+4)^2+(x2-6)^2+(x3+2)^2-4;
c2 = (x1+3)^2+(x2+4)^2+(x3+4)^2-4;
c3 = (x1-4)^2+(x2-0)^2+(x3-0)^2-5;
c4 = (x1+4)^2+(x2-2)^2+(x3-4)^2-5;
figure(1);clf;hold on;
% c1s = char(sdisplay(c1));
c2s = char(sdisplay(c2));
c3s = char(sdisplay(c3));
c4s = char(sdisplay(c4));
xran=[-dom dom -dom dom -dom dom];
% smrplot(c1s,0,xran,[400 100],'k-');
smrplot(c2s,0,xran,[400 100],'k-');
smrplot(c3s,0,xran,[400 100],'k-');
smrplot(c4s,0,xran,[400 100],'k-');
refreshdata; drawnow;
% % largest level set
c_0 = 0.1;
cc = 1;
iter = 1;
while abs(double(c_0)-double(cc)) >= 1e-6
    iter = iter + 1;
    if iter ~= 1
        c_0 = cc;
    end
    %% fix cc search for u, L1
    [uu1 u1c u1v] = polynomial([x1 x2 x3], 4, 0);
    [uu2 u2c u2v] = polynomial([x1 x2 x3], 4, 0);
    [uu3 u3c u3v] = polynomial([x1 x2 x3], 4, 0);
    [L1 L1c L1v] = polynomial([x1 x2 x3], 4, 0);
    %%
    Vdot = jacobian(v, x1)*(f(1)+uu1)+jacobian(v, x2)*(f(2)+uu2)+jacobian(v, x3)*(f(3)+uu3);
    %%
    F = [sos(L1), sos(-Vdot-L1*(cc-v))];
    [sol,vv,QQ] = solvesos(F,1,[],[L1c;u1c;u2c;u3c]);
    %%
    solL1 = value(L1c)'*L1v;
    solu1 = value(u1c)'*u1v;
    solu2 = value(u2c)'*u2v;
    solu3 = value(u3c)'*u3v;
    
    % fix solu,solL1 search for vc
    sdpvar vc
    [uu1 u1c u1v] = polynomial([x1 x2 x3], 4, 0);
    [uu2 u2c u2v] = polynomial([x1 x2 x3], 4, 0);
    [uu3 u3c u3v] = polynomial([x1 x2 x3], 4, 0);
%     [L2 L2c L2v] = polynomial([x1 x2 x3], 4, 0);
    [L3 L3c L3v] = polynomial([x1 x2 x3], 4, 0);
    [L4 L4c L4v] = polynomial([x1 x2 x3], 4, 0);
    [L5 L5c L5v] = polynomial([x1 x2 x3], 4, 0);
    %%
    Vdot = jacobian(v, x1)*(f(1)+uu1)+jacobian(v, x2)*(f(2)+uu2)+jacobian(v,x3)*(f(2)+uu3);
%     F = [sos(L2),sos(L3),sos(L4),sos(L5),sos(-Vdot-solL1*(vc-v)),sos(-(vc-v)+c1*L2),sos(-(vc-v)+c2*L3),sos(-(vc-v)+c3*L4),sos(-(vc-v)+c4*L5),vc>0];
    F = [sos(L3),sos(L4),sos(L5),sos(-Vdot-solL1*(vc-v)),sos(-(vc-v)+c2*L3),sos(-(vc-v)+c3*L4),sos(-(vc-v)+c4*L5),vc>0];
%%
%     [sol,vv,QQ] = solvesos(F,-vc,[],[u1c;u2c;u3c;L2c;L3c;L4c;L5c;vc]);
    [sol,vv,QQ] = solvesos(F,-vc,[],[u1c;u2c;u3c;L3c;L4c;L5c;vc]);
    %%
    solv = value(vc);
    cc = solv;
    % smr variables
    solh = value(vc)-v;
    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-8 8 -8 8 -8 8];
    if mod(real(iter),2) == 0
        smrplot(hs,0,xran,[300 50],'g-');
    else
        smrplot(hs,0,xran,[300 50],'r-');
    end
    axis(xran)
    refreshdata; drawnow;
end

%%
figure(1);hold on;
% find barrier certificates
% init h with CLF
solh = value(vc)-v;
hs = char(sdisplay(solh));
xran=[-8 8 -8 8 -8 8];
smrplot(hs,0,xran,[300 50],'b-');
refreshdata; drawnow;
Control = [];
Barrier = [];
%%
L_u = 4; L_us = 4; L_au = 6; h_degree = 4;
%%
k_us = 4; k_au = 6; B_degree = 4; k_u = 4; gamma = 0;
%%
% for ii = 1:20
iter = 1;
while 1
    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-dom dom -dom dom -dom dom];
    if mod(real(iter),2) == 0
        smrplot(hs,0,xran,[300 50],'k-');
    else
        smrplot(hs,0,xran,[300 50],'c-');
    end
    iter = iter + 1;
    sdpvar htol
    % fix h search for u, L1, L2
    [u1 u1c u1v] = polynomial([x1 x2 x3], L_u, 0);
    [u2 u2c u2v] = polynomial([x1 x2 x3], L_u, 0);
    [u3 u3c u3v] = polynomial([x1 x2 x3], L_u, 0);
    [L1 L1c L1v] = polynomial([x1 x2 x3], L_au, 0);
    [L2 L2c L2v] = polynomial([x1 x2 x3], L_au, 0);
    %%
    hdot = jacobian(solh, x1)*(f(1)+u1)+jacobian(solh, x2)*(f(2)+u2)+jacobian(solh, x3)*(f(3)+u3);
    Vdot = jacobian(v, x1)*(f(1)+u1)+jacobian(v, x2)*(f(2)+u2)+jacobian(v, x3)*(f(3)+u3);
    %%
    F = [sos(L1),sos(L2),sos(-Vdot+L1*(-solh)),sos(hdot+gamma*solh+L2*(-solh)-htol),htol>=0];
    [sol,vv,QQ] = solvesos(F,-htol,[],[L1c;L2c;u1c;u2c;u3c]);
    %%
    solL1 = value(L1c)'*L1v;
    solL2 = value(L2c)'*L2v;
    solu1 = value(u1c)'*u1v;
    solu2 = value(u2c)'*u2v;
    solu3 = value(u3c)'*u3v;
    Control = [Control; [solu1 solu2 solu3]];
    %%
    % fix u, L1, L2, search for h, L3, L4
    [h hc hv] = polynomial([x1 x2 x3], h_degree, 0);
    [L3 L3c L3v] = polynomial([x1 x2 x3], L_us, 0);
    [L4 L4c L4v] = polynomial([x1 x2 x3], L_us, 0);
    [L5 L5c L5v] = polynomial([x1 x2 x3], L_us, 0);
    %%
    hdot = jacobian(h, x1)*(f(1)+solu1)+jacobian(h, x2)*(f(2)+solu2)+jacobian(h, x3)*(f(3)+solu3);
    Vdot = jacobian(v, x1)*(f(1)+solu1)+jacobian(v, x2)*(f(2)+solu2)+jacobian(v, x3)*(f(3)+solu3);
    %%
    F = [sos(L3),sos(L4),sos(L5),sos(-Vdot+solL1*(-h)),sos(hdot+gamma*h+solL2*(-h)),sos(-h+c2*L3),sos(-h+c3*L4),sos(-h+c4*L5)];
    [sol,vv,QQ] = solvesos(F,-(hc(1)+hc(5)+hc(7)+hc(10)+hc(21)+hc(23)+hc(25)+hc(30)+hc(32)+hc(35)),[],[L3c;L4c;L5c;hc]);
    %%
    solh = value(hc)'*hv;
    B = solh;
    Barrier = [Barrier; solh];
    %%
    sdisplay(solh)
    figure(1);hold on;
    hs = char(sdisplay(solh));
    xran=[-dom dom -dom dom -dom dom];
    if mod(real(iter),2) == 0
        smrplot(hs,0,xran,[300 50],'m-');
    else
        smrplot(hs,0,xran,[300 50],'b-');
    end
    axis(xran)
    refreshdata; drawnow;
    
    %%
    sdpvar htol2
    [L1 L1c L1v] = polynomial([x1 x2 x3], k_au, 0);
    [L2 L2c L2v] = polynomial([x1 x2 x3], k_au, 0);
    [L3 L3c L3v] = polynomial([x1 x2 x3], k_us, 0);
    [L4 L4c L4v] = polynomial([x1 x2 x3], k_us, 0);
    [L5 L5c L5v] = polynomial([x1 x2 x3], k_us, 0);
    [h hc hv] = polynomial([x1 x2 x3], B_degree, 0);
    %%
    hdot = jacobian(h,x1)*(f(1)+solu1)+jacobian(h,x2)*(f(2)+solu2)+jacobian(h,x3)*(f(3)+solu3);
    %%
    F = [sos(L1),sos(L2),sos(L3),sos(L4),sos(L5),sos(h-L1*B-htol2),sos(hdot+L2*B+gamma*B-htol2),htol2>=0,sos(-h+c2*L3),sos(-h+c3*L4),sos(-h+c4*L5)];
    [sol,vv,QQ] = solvesos(F,-htol2,[],[L1c;L2c;L3c;L4c;L5c;hc;htol2]);
    %%
    solh = value(hc)'*hv;
    sdisplay(solh)
    % smr variables
    figure(real(100));clf;hold on;
    xran=[-dom dom -dom dom -dom dom];
    hs = char(sdisplay(solh));
    if mod(real(iter),2) == 0
        smrplot(hs,0,xran,[300 50],'r-');
    else
        smrplot(hs,0,xran,[300 50],'g-');
    end
    refreshdata; drawnow;
end

