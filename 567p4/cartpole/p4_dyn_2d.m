function dxdt = dyn_controller_paper_1d(~,x)
%DYN2D time invariant 2D nonlinear time continous dynamics
% In:
%    t      ~       time (not used)
%    x    2 x N     state
% Out:
%    dxdt 2 x N     state derivative
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

% dxdt(1,:) = x(2,:)-1.0*x(1,:);
% term = 0;
% deg = length(coe);
% for i=1:length(coe)
%     term = term+(x(1,:).^(deg-i+1)).*double(coe(1));
% end
% dxdt(2,:) = x(1,:)^2*x(2,:)+term;

g = 10;
l = 10;
m = 1;
M = 3;
dxdt(1,:) = x(1,:);
dxdt(2,:) = -(13*m*l*sin(x()));
dxdt(3,:) = x(3,:);
dxdt(4,:) = 

end