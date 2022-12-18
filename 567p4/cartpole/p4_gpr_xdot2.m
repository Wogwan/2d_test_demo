function [mean1,hyp1,delta,rmse] = p4_gpr_xdot2(x,y,xtest,ytest,it,noise,poly_deg,plott)
%%SOS_GP_CHECK_XDOT2 generates the xdot2 value [MEAN] to construct the polynomial function
% In:
%     Xtr_0      double   400  x  2   Input training X
%     dXtr_0     double   400  x  2   Input training Y
% Out:
%     mean1      double   6  x  1    Learned kernel with polynomial m
% Copyright (c) by Huang Hejun (CUHK) under BSD License
% Last modified: Huang Hejun 2021-05

tic
rng(1)
syms x1 x2;
dxdt2 = 0;
X = [1];
%% Mean function
m1 = {@meanPoly,poly_deg};
hyp_m1 = zeros([1*poly_deg 1]);
m2 = {@meanConst}; 
hyp_m2 = noise;
meanfunc = {'meanSum',{m2,m1}}; hyp.mean = [hyp_m2; hyp_m1];
%%
sf = 0.2; ell = 0.4;
cov1 = {@covSEiso}; hyp_cov1 = log([ell;sf/2]);
covfunc = cov1; hyp.cov = hyp_cov1;

%% Lik function
likfunc = @likGauss; sn = noise; hyp.lik = log(sn);

%% Inf function
inf_func = {'infGaussLik'};

%%
hyp1 = minimize(hyp, @gp, -it, inf_func, meanfunc, covfunc, likfunc, x, y);
%%
[ymu, ys2, fmu, fs2] = gp(hyp1, inf_func, meanfunc, covfunc, likfunc, x, y, xtest);

mean1 = hyp1.mean;
toc
%%
pvar x1 x2
X1 = p2s(monomials(x1,(1:poly_deg)));
for i = 1:length(X1)
    X = [X; X1(i)];
end
for i = 1:length(X)
    dxdt2 = dxdt2 + X(i)*mean1(i);
end
delta = 0;

%%
n = size(xtest,1);

%% Gathered
    error = abs(ytest-ymu);
    rmse = sqrt(sum(error.^2)/n);
    time1 = toc;
    %%
    syms x1 x2
    dXtr_3 = [];
    for num = 1:length(xtest(:,1))
        num_0 = vpa(subs(dxdt2,[x1],[xtest(num,1)]));
        dXtr_3 = [dXtr_3; num_0];
    end
    dXtr_3 = reshape(double(dXtr_3),1,[])';
    error_2 = abs(ytest-dXtr_3);
    rmse_poly = sqrt(sum(error_2.^2)/n);

if plott==1
    figure(801);clf;hold on;
    hold on;
    le1 = plot((1:n)',ymu,'k*','LineWidth',0.8);
    le2 = plot((1:n)',ytest, 'r+', 'LineWidth',0.7);
    le3 = plot((1:n)',dXtr_3,'bo','LineWidth',0.7);
    set(gca,'FontSize',24,'Fontname','Times');
    set(gca,'Box','on');
    ax = gca;
    ax.BoxStyle = 'full';
    ax.LineWidth = 1.2;
    xlabel('$x_1$','Interpreter','latex','Fontsize',24,'Fontname','Times');
    ylabel('$d_{\xi}$','Interpreter','latex','Fontsize',24,'Fontname','Times');
    legend([le2;le1;le3],{'Data','SE Kernel','Polynomial'},'Interpreter','latex','Orientation','horizon');
end
end