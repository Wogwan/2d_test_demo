function [SOLu,SOL1,SOL2] = sos_function_1(k,solh,V,mm,gamma)
    pvar x1 x2 u htol epsi;
    x = [x1;x2];
    % Create corresponding decision variable
%     [L1,L1_Q] = polydecvar('L1_w',monomials(x,0:k)); % L1 sos decision variables
%     [L2,L2_Q] = polydecvar('L2_w',monomials(x,0:k)); % L1 sos decision variables
    [L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:k)); % L1 sos decision variables
    [L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:k)); % L1 sos decision variables
    [u1,u1_Q] = polydecvar('u1_w',monomials(x,0:k)); % u1 sos decision variables
    
    if mm > 1
        [~,~]=pcontour(solh,0,domain,'r'); hold on;             % Plot the original Lyapunov sublevel set
        trace_Q1 = trace_Q;
        hdot = jacobian(solh,x1)*x2 + jacobian(solh, x2)*(-x1 + u1);
%         Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + u1);
        Vdot = jacobian(-solh, x1)*x2 + jacobian(-solh, x2)*(-x1 + u1);
    else
        hdot = jacobian(solh,x1)*x2 + jacobian(solh, x2)*(-x1 + u1);
        Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + u1);
    end

%%
    % Constrain :
    %     sosconstr(1) = -Vdot >= L1*solh;
    %     sosconstr(1) = -Vdot  >= L1*solh;
    %     sosconstr(2) = hdot+gamma*solh >=  L2*solh+htol;
    %     sosconstr(3) = htol >= 0;
    %     sosconstr(4) = L1 >= 0;
    %     sosconstr(5) = L2 >= 0;

    sosconstr(1) = L1 >= 0;
    sosconstr(2) = L2 >= 0;
    if mm == 1
        sosconstr(3) = -Vdot>= L1*solh;
        sosconstr(4) = hdot+gamma*solh >=  L2*solh+htol;
    else
        sosconstr(3) = -Vdot >= L1*solh;
        sosconstr(4) = hdot+gamma*solh >=  L2*solh+htol;    
    end


    %     sosconstr(2) = hdot+gamma*solh >=  L2*solh;
    %     sosconstr(3) = L1 >= 0;
    %     sosconstr(3) = L2 >= 0;

%% Set objection
    obj = -htol;

    % Solve feasibility problem
    opts = sosoptions;
    opts.form = 'kernel';
    opts.solver = 'mosek';
%     opts.solver = 'sedumi';
%     opts.solver = 'sdpam';
    [info,dopt] = sosopt(sosconstr,x,obj,opts);
    %     [info,dopt] = sosopt(sosconstr,x,opts);

    % Create output
    if info.feas
        SOL1 = subs(L1,dopt);
        SOL2 = subs(L2,dopt);
        SOLu = subs(u1,dopt);
    else
        kk = 0;
        fprintf('L can not find.====== ');
    end
end