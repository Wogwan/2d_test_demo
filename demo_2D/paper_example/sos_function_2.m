function [solh,trace_Q,kk]=sos_function_2(f,k,SOLu,SOL1,SOL2,gamma,mm,V,C,dom,gg,L_us)
kk = 1;
domain = [-dom dom -dom dom];
pvar x1 x2 htol epsi;
x = [x1;x2];
[h,hc] = polydecvar('h_w',monomials(x,0:k)); % L1 sos decision variables
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:L_us/2)); % L1 sos decision variables
[L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:L_us/2)); % L1 sos decision variables
[L5,L5_Q] = sosdecvar('L5_w',monomials(x,0:L_us/2)); % L1 sos decision variables
% [L6,L6_Q] = sosdecvar('L6_w',monomials(x,0:k/2)); % L1 sos decision variables
if mm == 0
    hdot = jacobian(h, x1)*f(1) + jacobian(h, x2)*(f(2)+gg*SOLu);
    Vdot = jacobian(V, x1)*f(1) + jacobian(V, x2)*(f(2)+gg*SOLu);
else
    hdot = jacobian(h, x1)*f(1) + jacobian(h, x2)*(f(2)+gg*SOLu);
    Vdot = jacobian(V, x1)*f(1) + jacobian(V, x2)*(f(2)+gg*SOLu);
end
%% Constraint:
if mm == 0
    pconstr_1 = -Vdot-SOL1*h >= 0;
    pconstr_2 = hdot+gamma*h-SOL2*h >= 0;
else
    pconstr_1 = -SOL1*h-Vdot >= 0;
    pconstr_2 = hdot+gamma*h-SOL2*h >= 0;
end
pconstr_3 = -h+C(1)*L3 >= 0;
pconstr_4 = -h+C(2)*L4 >= 0;
pconstr_5 = -h+C(3)*L5 >= 0;
pconstr_6 = L3 >= 0;
pconstr_7 = L4 >= 0;
pconstr_8 = L5 >= 0;

pconstr = [pconstr_6;pconstr_7;pconstr_8;pconstr_1;pconstr_2;pconstr_3;pconstr_4;pconstr_5];

%% Set objection
if k==2
    obj = -(hc(1)+hc(4)+hc(6));
elseif k==4
    obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15));
elseif k==6
    obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28));
elseif k==8
    obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45));
elseif k==10
    obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45)+hc(56)+hc(58)+hc(60)+hc(62)+hc(64)+hc(66));
elseif k==12
    obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45)+hc(56)+hc(58)+hc(60)+hc(62)+hc(64)+hc(66)+hc(79)+hc(81)+hc(83)+hc(85)+hc(87)+hc(89)+hc(91));
elseif k==14
    obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45)+hc(56)+hc(58)+hc(60)+hc(62)+hc(64)+hc(66)+hc(79)+hc(81)+hc(83)+hc(85)+hc(87)+hc(89)+hc(91)+hc(106)+hc(108)+hc(110)+hc(112)+hc(114)+hc(116)+hc(118)+hc(120));
elseif k==16
    obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45));
else
    fprintf('Pleaes enter suitable order of Barrier certificate.====== ');
end
%     obj = -sum(diag(hc));
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
figure(12);hold on;
%% Create output
if info.feas
    solh = subs(h,dopt);
    Q = subs(-obj,dopt);
    trace_Q = Q;
    [~,~]=pcontour(solh,0,domain,'g'); hold on;             % Plot the original Lyapunov sublevel set
    refreshdata; drawnow;
else
    kk = 0;
    solh = 0;
    trace_Q = 0;
    Q = 0;
    fprintf('Barrier Certificate can not find.======\n');
    return;
end
end