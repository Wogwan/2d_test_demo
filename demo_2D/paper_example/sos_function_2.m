function [solh,trace_Q,kk]=sos_function_2(iter,f,k,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us)
kk = 1;
domain = [-dom dom -dom dom];
pvar x1 x2 htol epsi;
x = [x1;x2];
%%
[h,hc] = polydecvar('h_w',monomials(x,0:k)); % L1 sos decision variables
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:L_us/2)); % L1 sos decision variables
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:L_us/2)); % L1 sos decision variables
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:L_us/2)); % L1 sos decision variables
%%
hdot = jacobian(h, x1)*(f(1)+gg(1)*SOLu1)+jacobian(h, x2)*(f(2)+gg(2)*SOLu2);
Vdot = jacobian(V, x1)*(f(1)+gg(1)*SOLu1)+jacobian(V, x2)*(f(2)+gg(2)*SOLu2);
%% Constraint:
pconstr_6 = L1 >= 0;
pconstr_7 = L2 >= 0;
pconstr_8 = L3 >= 0;
pconstr_1 = -Vdot-SOL1*h >= 0;
pconstr_2 = hdot+gamma*h-SOL2*h >= 0;
pconstr_3 = -h+C(1)*L1 >= 0;
pconstr_4 = -h+C(2)*L2 >= 0;
pconstr_5 = -h+C(3)*L3 >= 0;
% pconstr = [pconstr_6;pconstr_7;pconstr_8;pconstr_3;pconstr_4;pconstr_5;pconstr_1;pconstr_2];
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
else
    fprintf('Pleaes enter suitable order of Barrier certificate.====== ');
end
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,obj,opts);
figure(12);hold on;
%% Create output
if info.feas
    kk = 1;
    solh = subs(h,dopt);
    trace_Q = subs(-obj,dopt);
    k = ['r','g','b','m','c','k','y'];
    if mod(iter,7) == 0
        [~,~]=pcontour(solh,0,domain,k(7)); hold on;             % Plot the original Lyapunov sublevel set
    else
        [~,~]=pcontour(solh,0,domain,k(mod(iter,7))); hold on;             % Plot the original Lyapunov sublevel set
    end
    refreshdata; drawnow;
else
    kk = 0;
    solh = 0;
    trace_Q = 0;
    fprintf('Barrier Certificate can not find.======\n');
    return;
end
end