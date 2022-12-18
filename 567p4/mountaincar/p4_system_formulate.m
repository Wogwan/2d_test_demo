function sys = p4_system_formulate

%% System hyperparameters
sys.dom = 20; 
sys.plott = 0; % Plot figure
sys.domain = [-sys.dom sys.dom -sys.dom sys.dom];
sys.g = 10;
sys.l = 10;
sys.m = 1;
% sys.gg = [1;1/(sys.m*(sys.l^2))];
sys.gg = [0;1];
%% The first dimension
syms x1 x2
sys.f2_npoly = -sys.g/sys.l*sin(x1);
%%
pvar x1 x2
sys.V = [
    x1^2+x2^2+x1*x2
    x1^2+x2^2+x1*x2+x1^4+x2^4
    1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2
    1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1
    ];
% x1 is using the standard 1 rad
C1 = x1 + 3;
C2 = -x1 + 3;
% C1 = (x1+4)^2+(x2+4)^2-3;
% C2 = (x1-4)^2+(x2-4)^2-3;
% sys.us_region = [C1;C2];
C3 = x2 + 10;
C4 = -x2 + 10;
sys.us_region = [C1;C2;C3;C4];
end