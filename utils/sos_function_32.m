function solu = sos_function_32(k,h,gamma,mm,V,C,V0)
    
    record = 0;
    domain = [-8 8 -8 8];   
    pvar x1 x2;   
    x = [x1;x2];
    [u,hc] = polydecvar('u_w',monomials(x,0:k)); % L1 sos decision variables

%     [L7,L7_Q] = sosdecvar('L7_w',monomials(x,0:k)); % L1 sos decision variables
%     [L8,L8_Q] = sosdecvar('L8_w',monomials(x,0:k)); % L1 sos decision variables
%     [L9,L9_Q] = sosdecvar('L9_w',monomials(x,0:k)); % L1 sos decision variables
%     [L10,L10_Q] = sosdecvar('L10_w',monomials(x,0:k)); % L1 sos decision variables

    hdot = jacobian(h,x1)*x2 + jacobian(h, x2)*(-x1 + u);
    Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + u);
    
%     if mm > 1
% % Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + u1);
%         Vdot = jacobian(solh, x1)*x2 + jacobian(solh, x2)*(-x1 + SOLu);
%     else
%         Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + SOLu);
%     end
    
%% Constraint:

    pconstr(1) = h >= x1^2+x2^2;
    pconstr(2) = hdot+gamma*h >= 0;
    pconstr(3) = -Vdot >= x1^2+x2^2;

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
        solu = subs(u,dopt);
        record = solu
        U = subs(obj,dopt);
        trace_U = U
        [~,~]=pcontour(solh,0,domain,'g'); hold on;             % Plot the original Lyapunov sublevel set
        refreshdata; drawnow;
    else
        kk = 0;
        solu = record;
        trace_U = record;
        fprintf('Controller is found.======\n');
        return;
    end
end