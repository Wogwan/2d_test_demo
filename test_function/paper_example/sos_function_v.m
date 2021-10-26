function [solu,solL,kk]=sos_function_v(f,k,V,cc,dom)
    kk = 1;
    domain = [-dom dom -dom dom];   
    pvar x1 x2 u;   
    x = [x1;x2];
%     [h,hc] = sosdecvar('h_w',monomials(x,0:k/2)); % L1 sos decision variables
    [u,uc] = polydecvar('u_w',monomials(x,0:k)); % L1 sos decision variables
    
    [L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k/2)); % L1 sos decision variables
%     [L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:k/2)); % L1 sos decision variables
%     [L5,L5_Q] = sosdecvar('L5_w',monomials(x,0:k/2)); % L1 sos decision variables
%     [L6,L6_Q] = sosdecvar('L6_w',monomials(x,0:k/2)); % L1 sos decision variables

    Vdot = jacobian(V, x1)*(0.1*x1^2+1*x2)+ jacobian(V, x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*u);
%     hdot = jacobian(h,x1)*x2 + jacobian(h, x2)*(-x1 + SOLu);
%     Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + SOLu);    
%     if mm > 1
% %         Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + u1);
%         Vdot = jacobian(V, x1)*(0.1*x1^2+1*x2) + jacobian(V, x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*SOLu);
%     else
%         Vdot = jacobian(V, x1)*(0.1*x1^2+1*x2) + jacobian(V, x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*SOLu);
%     end
    
%% Constraint:
        pconstr_1 = L3 >= 0;
        pconstr_2 = -Vdot-L3*(cc-V) >= 0;
        pconstr = [pconstr_1; pconstr_2];
        
%% Set objection
%     obj = -sum(diag(hc));

%% Solve feasibility problem
    opts = sosoptions;
    opts.form = 'kernel';
    opts.solver = 'mosek';
%     opts.solver = 'sedumi';
%     opts.solver = 'sdpam';
    %     [info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);

    [info,dopt] = sosopt(pconstr,x,opts);

    % Create output
    if info.feas
        solu = subs(u,dopt);
        solL = subs(L3,dopt);
    else
        kk = 0;
        solu  = 0;
        solL = 0;
        fprintf('Barrier Certificate can not find.======\n');
        return;
    end
end