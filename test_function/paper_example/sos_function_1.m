function [SOLu,SOL1,SOL2, kk] = sos_function_1(f,k,solh,V,mm,gamma,dom)
pvar x1 x2 u htol epsi;
x = [x1;x2];
% Create corresponding decision variable
%     [L1,L1_Q] = polydecvar('L1_w',monomials(x,0:k)); % L1 sos decision variables
%     [L2,L2_Q] = polydecvar('L2_w',monomials(x,0:k)); % L1 sos decision variables
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:k/2)); % L1 sos decision variables
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:k/2)); % L1 sos decision variables
[u,u1_Q] = polydecvar('u1_w',monomials(x,0:k)); % u1 sos decision variables

figure(12); hold on; domain = [-dom dom -dom dom];
%     if mm > 1
%         [~,~]=pcontour(solh,0,domain,'r'); hold on;             % Plot the original Lyapunov sublevel set
%         hdot = jacobian(solh,x1)*f(1) + jacobian(solh, x2)*f(2);
% %         Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + u1);
%         Vdot = jacobian(-solh, x1)*f(1) + jacobian(-solh, x2)*f(2);
%     else
%         hdot = jacobian(solh, x1)*(0.1*x1^2+1*x2) + jacobian(solh, x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*u);
%         Vdot = jacobian(V, x1)*(0.1*x1^2+1*x2) + jacobian(V, x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*u);
%     end
hdot = jacobian(solh, x1)*f(1) + jacobian(solh, x2)*(f(2)+gg*u);
Vdot = jacobian(V, x1)*(0.1*x1^2+1*x2) + jacobian(V, x2)*(0.1*x1*x2-0.2*x1+(1+x1^2)*u);

%%
% Constrain :
%     sosconstr(1) = -Vdot >= L1*solh;
%     sosconstr(1) = -Vdot  >= L1*solh;
%     sosconstr(2) = hdot+gamma*solh >=  L2*solh+htol;
%     sosconstr(3) = htol >= 0;
%     sosconstr(4) = L1 >= 0;
%     sosconstr(5) = L2 >= 0;

sosconstr_1 = L1 >= 0;
sosconstr_2 = L2 >= 0;
if mm == 1
    sosconstr_3 = -Vdot>= L1*solh;
    sosconstr_4 = hdot+gamma*solh >=  L2*solh+htol;
else
    sosconstr_3 = -Vdot >= L1*solh;
    sosconstr_4 = hdot+gamma*solh >=  L2*solh+htol;
end
sosconstr_5 = htol >=0;

sosconstr = [sosconstr_1;sosconstr_2;sosconstr_3;sosconstr_4;sosconstr_5];
%     sosconstr = [sosconstr_3;sosconstr_4;sosconstr_5;sosconstr_1;sosconstr_2];

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
    SOLu = subs(u,dopt);
    kk = 1;
else
    kk = 0;
    fprintf('L can not find.====== ');
end
end