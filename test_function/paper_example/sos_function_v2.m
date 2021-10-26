function [cc,kk]=sos_function_v2(f,k,V,C,dom,solL)
    kk = 1;
    domain = [-dom dom -dom dom];   
    pvar x1 x2 cc;   
    x = [x1;x2];

%     [h,hc] = sosdecvar('h_w',monomials(x,0:k/2)); % L1 sos decision variables
    
    [u,uc] = polydecvar('u_w',monomials(x,0:k)); % L1 sos decision variables
    [L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:k/2)); % L1 sos decision variables


    Vdot = jacobian(V, x1)*(0.1*x1^2+1*x2) + jacobian(V, x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*u);
 
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
        pconstr_2 = L4 >= 0;
        pconstr_1 = -Vdot-solL*(cc-V) >= 0;
        pconstr_3 = -(cc-V)+C(1)*L4 >= 0;
        pconstr_4 = cc >= 0;
    
        pconstr = [pconstr_2;pconstr_1;pconstr_3;pconstr_4];              

%% Set objection
    obj = -cc;

%% Solve feasibility problem
    opts = sosoptions;
    opts.form = 'kernel';
    opts.solver = 'mosek';
    [info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
    % [info,dopt] = sosopt(pconstr,x,obj);
    
    figure(11);hold on;
    % Create output
    if info.feas
        cc = subs(cc,dopt);
        [~,~]=pcontour(V,double(cc),domain,'g'); hold on;             % Plot the original Lyapunov sublevel set
        refreshdata; drawnow;
    else
        kk = 0;
        cc = 0;
        fprintf('Barrier Certificate can not find.======\n');
        return;
    end
end