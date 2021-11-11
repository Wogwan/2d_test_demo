function [solh,trace_Q,kk]=sos_function_2_3D(iter,f,k,SOLu1,SOLu2,SOLu3,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id)
kk = 1;
domain = [-dom dom -dom dom -dom dom];
pvar x1 x2 x3;
x = [x1;x2;x3];
[h,hc] = polydecvar('h_w',monomials(x,0:k)); 
[L3,L3_Q] = polydecvar('L3_w',monomials(x,0:L_us));
[L4,L4_Q] = polydecvar('L4_w',monomials(x,0:L_us)); 
[L5,L5_Q] = polydecvar('L5_w',monomials(x,0:L_us)); 
[L6,L6_Q] = polydecvar('L6_w',monomials(x,0:L_us)); 
hdot = jacobian(h, x1)*(f(1)+gg(1)*SOLu1)+jacobian(h, x2)*(f(2)+gg(2)*SOLu2)+jacobian(h, x3)*(f(3)+gg(3)*SOLu3);
Vdot = jacobian(V, x1)*(f(1)+gg(1)*SOLu1)+jacobian(V, x2)*(f(2)+gg(2)*SOLu2)+jacobian(V, x3)*(f(3)+gg(3)*SOLu3);
%% Constraint:
pconstr_1 = -Vdot-SOL1*h >= 0;
pconstr_2 = hdot+gamma*h-SOL2*h >= 0;
pconstr_41 = L3 >= 0;
pconstr_42 = L4 >= 0;
pconstr_43 = L5 >= 0;
pconstr_44 = L6 >= 0;
pconstr_31 = -h+C(1)*L3 >= 0;
pconstr_32 = -h+C(2)*L4 >= 0;
pconstr_33 = -h+C(3)*L5 >= 0;
pconstr_34 = -h+C(4)*L6 >= 0;
pconstr = [pconstr_41;pconstr_42;pconstr_43;pconstr_44;pconstr_1;pconstr_2;pconstr_31;pconstr_32;pconstr_33;pconstr_34];

%% Set objection
if k==2
    obj = -(hc(1)+hc(5)+hc(7)+hc(10));
elseif k==4
    obj = -(hc(1)+hc(5)+hc(7)+hc(10)+hc(21)+hc(23)+hc(25)+hc(30)+hc(32)+hc(35));
elseif k==6
    obj = -(hc(1)+hc(5)+hc(7)+hc(10)+hc(21)+hc(23)+hc(25)+hc(30)+hc(32)+hc(35)+hc(57)+hc(59)+hc(61)+hc(63)+hc(70)+hc(72)+hc(74)+hc(79)+hc(81)+hc(84));
elseif k==8
    obj = -(hc(1)+hc(5)+hc(7)+hc(10)+hc(21)+hc(23)+hc(25)+hc(30)+hc(32)+hc(35)+hc(57)+hc(59)+hc(61)+hc(63)+hc(70)+hc(72)+hc(74)+hc(79)+hc(81)+hc(84)+hc(121)+hc(123)+hc(125)+hc(127)+hc(129)+hc(131)+hc(138)+hc(140)+hc(142)+hc(144)+hc(151)+hc(153)+hc(155)+hc(160)+hc(162)+hc(165));
else
    fprintf('Pleaes enter suitable order of Barrier certificate.====== ');
end
%     obj = -sum(diag(hc));
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,obj,opts);
figure(figure_id);hold on;
%% Create output
if info.feas
    solh = subs(h,dopt);
    Q = subs(-obj,dopt);
    trace_Q = Q;
    inH = patch(pcontour3(solh,0,domain,'k'));
    if mod(iter,2) == 1 
        set(inH, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'm','LineStyle','-','LineWidth',0.7 ); hold on;
        view(-30,20);hold on;
    else
        set(inH, 'EdgeAlpha',0.1,'FaceColor', 'none', 'EdgeColor', 'b','LineStyle','-','LineWidth',0.7 ); hold on;    
        view(-30,20);hold on;
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