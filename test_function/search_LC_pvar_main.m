function [L_t, k_mark] = search_LC_pvar_main(f,V,C,k)
%%SEARCH_L_PVAR_MAIN generates the L1 and L2 factor based on Algorithms step2 with SOSOPT toolbox
% In:
%     f     pvar-SOS    2  x  1   Input tracking system
%     V     pvar-SOS    1  x  1   Input existed Lyapunov function
%     h     pvar-SOS    1  x  1   Input existed Barrier certificate
% Out:
%     s1_sol     pvar-SOS   1  x  1    Factor L1
%     s2_sol     pvar-SOS   1  x  1    Factor L2
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

    pvar x1 x2;
    x = [x1;x2];
    k_mark = 0;
%%%%%%%%%%%%%%%%%%%%%%%
%     f = [x2;(1-x1^2)*x2-x1];
%     V = x1^2+x2^2+x1*x2;
%     C = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%
    dVdt = jacobian(V,x)*f;
    % Create corresponding decision variable
    [L1,L1_Q] = polydecvar('L1_w',monomials(x,0:k)); % L1 sos decision variables

    % Constraint 1: -Vdot-L*(c-V)*h  in SOS
    pconstr = -dVdt >= L1*(C-V);

    % Solve feasibility problem
    opts = sosoptions;
    opts.form = 'kernel';
    opts.solver = 'mosek';
    % opts.solver = 'sedumi';
    % opts.solver = 'sdpam';
    [info,dopt] = sosopt(pconstr,x,opts);
    % [info,dopt] = sosopt(pconstr,x,obj);

    % Create output
    if info.feas
        L_t = subs(L1,dopt);
    else
        L_t = [];
        fprintf('L can not find.====== ');
        k_mark = 1;
    end
end