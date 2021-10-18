pvar x1 x2 u;
x = [x1;x2];
k = 2;
%%%%%%%%%%%%%%%%%%%%%%%
% f = [0.89656063674407948660416423081188*x1^3-1.1066459905781698830340431527475*x1^2-0.94292523209241957404813661014487*x1+0.0040335440741594277488935027520256*x2^3-0.012631023336381969057740093376196*x2^2+0.018814164112983571691684048232673*x2-0.74128820430649877692985683097504*x3^3+0.57353888754769488667051291486132*x3^2-0.15011920973339787366285236203112*x3+0.026795052581369493971408246579813; 
%     -1.0*x2*x1^3-1.0*x2; 
%     -0.28781069452828345056616399233462*x1^4+0.44794385930119814953620505093568*x1^3-1.0*x1^2*x3-0.029564576378269480372296129644383*x1^2-0.13939107187891514039179696737847*x1-0.00522377410817732781150857235275*x2^3+0.00072167832110234188432162927284708*x2^2+0.058049930474332787910807240905342*x2+0.56325444655680989569646044401452*x3^3-0.088835620688520258725340283945116*x3^2-0.4943791483103421868783300396899*x3-0.0047284451719816766868120794242714];
%%%%%%%%%%%%%%%%%%%%%%%%
% [u,u_Q] = polydecvar('u_w',monomials(x,0:k)); % L1 sos decision variables
[L1,L1_Q] = polydecvar('L1_w',monomials(x,0:k)); % L1 sos decision variables
[L2,L2_Q] = polydecvar('L2_w',monomials(x,0:k)); % L1 sos decision variables
[V,V_Q] = polydecvar('V',monomials(x,0:k)); % L1 sos decision variables

u = x1 + x2;
f = [x2; -x1 + u];
Vdot = jacobian(V, x1)*x2 + jacobian(V, x2)*(-x1 + u);

dVdt = jacobian(V,x)*f;
% Constraint 1: -Vdot-L*(c-V)*h  in SOS
sosconstr(1) = L1 >= 0;
sosconstr(2) = L2 >= 0;
sosconstr(3) = V >= L1*(x1^2+x2^2);
sosconstr(4) = -Vdot >= L2*(x1^2+x2^2);

% Set objection
% obj = -C;

% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
% opts.solver = 'sedumi';
% opts.solver = 'sdpam';
% [info,dopt] = sosopt(pconstr,x,obj,opts);
[info,dopt] = sosopt(sosconstr,x,opts);

% Create output
if info.feas
    V_t = subs(V,dopt)
else
    V_t = [];
end