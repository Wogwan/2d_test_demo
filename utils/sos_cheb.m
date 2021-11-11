function [f2, f3] = sos_cheb(deg,sz)
%%SOS_CHEB generates the function from Chebyshev approximation in [g]
% In:
%     deg    double   1  x  1   Input approximation truncated series degree
%     sz     double   1  x  1   Input approximation in the given square region
% Out:
%      f2     syms     1 x 1     Polynomial function in systems
%      f3     syms     1 x 1     Polynomial function in systems

% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

% See corresponding [g] in SOS_MODEL.m.
% We can plot the approximation result with code below.

%%
%     clear;clf
%     deg = 4;                                                  % Target chebyshev truncated series degree
%     sz = 2;                                                   % Target chebyshev truncated series squared region
%%

    a = -sz; b = sz;
%     y2 = chebfun('-sin(x)+sin(x)^3',[a,b],'splitting','on'); % Modified here with different non-polynomial g
    y2 = chebfun('-cos(x)^2*sin(x)',[a,b],'splitting','on'); % Modified here with different non-polynomial g
    y_2 = minimax(y2,deg); c_2 = chebcoeffs(y_2);
    
%     y3 = chebfun('1-sqrt(sqrt((exp(x)*cos(x))^2))',[a,b],'splitting','on'); % Modified here with different non-polynomial g
    y3 = chebfun('1- (exp(2*x3)*cos(x3)^2)^(1/4)',[a,b],'splitting','on'); % Modified here with different non-polynomial g
    y_3 = minimax(y3,deg); c_3 = chebcoeffs(y_3);
        
%     Plot the chebyshev approximation result with original function
    figure(999);hold on;
    h1_1 = plot(y2,'b'); 
    h1_2 = plot(y_2,'k--');
    h2_1 = plot(y3,'r');
    h2_2 = plot(y_3,'g--');
    xlim([a b]);ylim([-3.5 3.5]);hold on;
    legend([h1_1,h1_2,h2_1,h2_2],{'$h1 = 1-\exp{(x^2)}$','$h1_2 = 1-\exp{(x^2)}$','$h2 = 1-\vert  \exp{(x)}\cos{(x)}\vert$','$h2_2 = 1-\vert \exp{(x)}\cos{(x)}\vert$'},'Interpreter','latex','location','northwest')
  
    syms x1 x3;
    T = chebyshevT([0:deg],x1);
    x2 = x1/sz;
    
    if length(c_2)~=(deg+1)
        c_2 = [c_2;0];
    end
    
    f2 = vpa(T*c_2);
    f2 = subs(f2,x1,x2);
    
    T = chebyshevT([0:deg],x3);
    f3 = vpa(T*c_3);
    x4 = x3/sz;
    f3 = subs(f3,x3,x4);
end

