function [solh, trace_Q, Q, kk]=sos_function_2(k,SOLu,SOL1,SOL2,gamma,mm,V,C)
    kk = 1;
    domain = [-8 8 -8 8];   
    pvar x1 x2 u htol epsi;   
    x = [x1;x2];
    [h,hc] = sosdecvar('h_w',monomials(x,0:k/2)); % L1 sos decision variables
%         [L3,L3_Q] = polydecvar('L3_w',monomials(x,0:k)); % L1 sos decision variables
%         [L4,L4_Q] = polydecvar('L4_w',monomials(x,0:k)); % L1 sos decision variables
%         [L5,L5_Q] = polydecvar('L5_w',monomials(x,0:k)); % L1 sos decision variables
%         [L6,L6_Q] = polydecvar('L6_w',monomials(x,0:k)); % L1 sos decision variables
    [L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k)); % L1 sos decision variables
    [L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:k)); % L1 sos decision variables
    [L5,L5_Q] = sosdecvar('L5_w',monomials(x,0:k)); % L1 sos decision variables
    [L6,L6_Q] = sosdecvar('L6_w',monomials(x,0:k)); % L1 sos decision variables

    hdot = jacobian(h,x1)*x2 + jacobian(h, x2)*(-x1 + SOLu);
%     Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + SOLu);    
    if mm > 1
%         Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + u1);
        Vdot = jacobian(solh, x1)*x2 + jacobian(solh, x2)*(-x1 + SOLu);
    else
        Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + SOLu);
    end
    
%% Constraint:
    %     pconstr(1) = -Vdot-SOL1*h >= 0;
    %     pconstr(2) = hdot+gamma*h-SOL2*h >= 0;
    %     pconstr(3) = h <= C1*L3;
    %     pconstr(4) = h <= C2*L4;
    %     pconstr(5) = h <= C3*L5;
    %     pconstr(6) = h <= C4*L6;
    %     pconstr(7) = L3 >= 0;
    %     pconstr(8) = L4 >= 0;
    %     pconstr(9) = L5 >= 0;
    %     pconstr(10) = L6 >= 0;

%%
        if mm == 1
            pconstr(1) = -Vdot-SOL1*h >= 0;
            pconstr(2) = hdot+gamma*h-SOL2*h >= 0;
        else
            pconstr(1) = SOL1*h-Vdot >= 0;
            pconstr(2) = hdot+gamma*h-SOL2*h >= 0;
        end
        pconstr(3) = h <= C(1)*L3;
        pconstr(4) = h <= C(2)*L4;
%         pconstr(5) = h <= C(3)*L5;
        pconstr(5) = L3 >= 0;
        pconstr(6) = L4 >= 0;
%         pconstr(8) = L5 >= 0;

%%
    %     pconstr(1) = L3 >= 0;
    %     pconstr(2) = L4 >= 0;
    %     pconstr(3) = L5 >= 0;
    %     pconstr(4) = -Vdot-SOL1*h >= 0;
    %     pconstr(5) = hdot+gamma*h-SOL2*h >= 0;
    %     pconstr(6) = h <= C1*L3;
    %     pconstr(7) = h <= C2*L4;
    %     pconstr(8) = h <= C3*L5;

%%
%     pconstr(1) = L3 >= 0;
%     pconstr(2) = L4 >= 0;
%     pconstr(3) = L5 >= 0;
%     pconstr(4) = L6 >= 0;
%     pconstr(5) = h <= C1*L3;
%     pconstr(6) = h <= C2*L4;
%     pconstr(7) = h <= C3*L5;
%     pconstr(8) = h <= C4*L6;
%     if mm == 1
%         pconstr(9) = -Vdot-SOL1*h >= 0;
%         pconstr(10) = hdot+gamma*h-SOL2*h >= 0;
%     else
%         pconstr(9) = SOL1*h-Vdot >= 0;
%         pconstr(10) = hdot+gamma*h-SOL2*h >= 0;
%     end

%%  Include C4
%     pconstr(1) = L3 >= 0;
%     pconstr(2) = L4 >= 0;
%     pconstr(3) = L5 >= 0;
%     pconstr(4) = h <= C1*L3;
%     pconstr(5) = h <= C2*L4;
%     pconstr(6) = h <= C3*L5;
%     if mm == 1
%         pconstr(7) = -Vdot-SOL1*h >= 0;
%         pconstr(8) = hdot+gamma*h-SOL2*h >= 0;
%     else
%         pconstr(7) = SOL1*h-Vdot >= 0;
%         pconstr(8) = hdot+gamma*h-SOL2*h >= 0;
%     end

%% Set objection
%     if k==2
%         obj = -(hc(1)+hc(4)+hc(6));
%     elseif k==4
%         obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15));
%     elseif k==6
%         obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28));
%     elseif k==8
%         obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45));
%     elseif k==10
%         obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45)+hc(56)+hc(58)+hc(60)+hc(62)+hc(64)+hc(66));
%     elseif k==12
%         obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45)+hc(56)+hc(58)+hc(60)+hc(62)+hc(64)+hc(66)+hc(79)+hc(81)+hc(83)+hc(85)+hc(87)+hc(89)+hc(91));
%     elseif k==14
%         obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45)+hc(56)+hc(58)+hc(60)+hc(62)+hc(64)+hc(66)+hc(79)+hc(81)+hc(83)+hc(85)+hc(87)+hc(89)+hc(91)+hc(106)+hc(108)+hc(110)+hc(112)+hc(114)+hc(116)+hc(118)+hc(120));
%     elseif k==16
%         obj = -(hc(1)+hc(4)+hc(6)+hc(11)+hc(13)+hc(15)+hc(22)+hc(24)+hc(26)+hc(28)+hc(37)+hc(39)+hc(41)+hc(43)+hc(45));  
%     else
%         fprintf('Pleaes enter suitable order of Barrier certificate.====== ');
%     end
    obj = -sum(diag(hc));

%% Solve feasibility problem
    opts = sosoptions;
    opts.form = 'kernel';
    opts.solver = 'mosek';
%     opts.solver = 'sedumi';
%     opts.solver = 'sdpam';
    %     [info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
    [info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
    % [info,dopt] = sosopt(pconstr,x,obj);
    
    figure(11);hold on;
    % Create output
    if info.feas
        solh = subs(h,dopt);
        record = solh
        Q = subs(-obj,dopt);
        trace_Q = Q
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