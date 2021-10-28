clear;clc;

sdpvar x1 x2 u
f = [x2-x1
    0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608
    ];
gu = x1+x2;
L1_factor = 4;

v = 1*x1^4+x1^2*x2^2+1*x2^4+1*x2^2+1*x2^2+x1*x2;
e = [0;0];
cc = 0.5;
gamma = 1;

% init h with CLF
solh = cc - v;
% safety constraints
c1 = (x1-1)^2+(x2-3)^2-1;
c2 = (x1+1)^2+(x2-3)^2-1;
c3 = (x1-3)^2+(x2-3)^2-1;
c4 = (x1+3)^2+(x2-3)^2-1;
c5 = (x1+2)^2+(x2+4)^2-1;
    figure(1);clf;hold on;
    c1s = char(sdisplay(c1));
    c2s = char(sdisplay(c2));
    c3s = char(sdisplay(c3));
    c4s = char(sdisplay(c4));
    c5s = char(sdisplay(c5));
    hinit = char(sdisplay(solh));
    xran=[-8 8 -8 8];
    smrplot(c1s,0,xran,[300 50],'b-');
    smrplot(c2s,0,xran,[300 50],'r-');
    smrplot(c3s,0,xran,[300 50],'r-');
    smrplot(c4s,0,xran,[300 50],'r-');
    smrplot(c5s,0,xran,[300 50],'r-');
    smrplot(hinit,0,xran,[300 50],'c-');
    axis(xran)

for ii = 1:10
    
    % fix cc search for u, L1
    [u uc uv] = polynomial([x1 x2], L1_factor, 0);
    [L1 L1c L1v] = polynomial([x1 x2], L1_factor, 0);
    
    % Vdot = jacobian(v, x1)*(-x2-3/2*x1^2-1.2*x1^3) + jacobian(v, x2)*(-x1 - u);
    Vdot = jacobian(v, x1)*f(1) + jacobian(v, x2)*(f(2)+gu*u);
    F = [sos(L1), sos(-Vdot - L1*(cc - v))];
    [sol,vv,QQ] = solvesos(F,1,[],[L1c;uc]);
    solL1 = value(L1c)'*L1v;
    solu = value(uc)'*uv;
    sdisplay(solL1)
    sdisplay(solu)
    
    % fix solu,solL1 search for vc
    sdpvar vc
    [u uc uv] = polynomial([x1 x2], L1_factor, 0);
    [L2 L2c L2v] = polynomial([x1 x2], L1_factor, 0);
    [L3 L3c L3v] = polynomial([x1 x2], L1_factor, 0);
    [L4 L4c L4v] = polynomial([x1 x2], L1_factor, 0);
    [L5 L5c L5v] = polynomial([x1 x2], L1_factor, 0);
    [L6 L6c L6v] = polynomial([x1 x2], L1_factor, 0);
    
    % Vdot = jacobian(v, x1)*(-x2-3/2*x1^2-1.2*x1^3) + jacobian(v, x2)*(-x1 - u);
    Vdot = jacobian(v, x1)*f(1) + jacobian(v, x2)*(f(2)+gu*u);

    %%
    F = [sos(L2), sos(-Vdot - solL1*(vc - v)), sos(-(vc-v)+c1*L2), vc>0];
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
