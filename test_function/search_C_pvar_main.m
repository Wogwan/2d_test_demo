function C_t = search_C_pvar_main(f,V,L)
%%SEARCH_C_PVAR_MAIN generates the sublevelset C based on Algorithms step1 with SOSOPT toolbox
% In:
%     f     pvar-SOS    2  x  1   Input tracking system
%     V     pvar-SOS    1  x  1   Input existed Lyapunov function
%     L     pvar-SOS    1  x  1   Polynomial Factor L
% Out:
%     C     double      1  x  1   Sub-level set C
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05


    pvar x1 x2 C;
    x = [x1;x2];
%%%%%%%%%%%%%%%%%%%%%%%
%     f = [x2-1.0*x1-0.000068554745018267748690732332761399;
%         -0.28403759803221195756606221038965*x1^4+0.52869001210367462074172766603977*x1^3+x1^2*x2+0.58373855464282150451861228056368*x1^2-0.82941162572600601876615655783098*x1*x2-1.080225579322203214042854831638*x1-0.0034675201912947187579683294700317*x2^2+0.21466335627848652057648135168296*x2+0.011335592411184153094405296854014];
%     V = 6*x1^2+6*x2^2-8*x1*x2;
%%%%%%%%%%%%%%%%%%%%%%%%
    dVdt = jacobian(V,x)*f;
    % Constraint 1: -Vdot-L*(c-V)*h  in SOS
    pconstr = -dVdt >= L*(C-V);
    % Constraint 2: Hdot+ga*h-S2*h-gam in SOS
    % pconstr(2) = dHdt+h>=L2*h;

    % Set objection
    obj = -C;

    % Solve feasibility problem
    opts = sosoptions;
    opts.form = 'kernel';
    opts.solver = 'mosek';
    % opts.solver = 'sedumi';
    % opts.solver = 'sdpam';
    [info,dopt] = sosopt(pconstr,x,obj,opts);
    % [info,dopt] = sosopt(pconstr,x,obj);

    % Create output
    if info.feas
        C_t = subs(C,dopt);
    else
        C_t = [];
    end
end