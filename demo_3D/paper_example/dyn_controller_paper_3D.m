function dxdt = dyn_controller_paper_3D(~,x)
%DYN2D time invariant 2D nonlinear time continous dynamics
% In:
%    t      ~       time (not used)
%    x    2 x N     state
% Out:
%    dxdt 2 x N     state derivative
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

% From main_total.m 25 line: dotx1 chebyshev approximation result with 4th order remainder
% dxdt(1,:) = x(2,:)-1.0*x(1,:);
% dxdt(2,:) = x(1,:)^2*x(2,:) - 0.38410921204577513909939057157317*x(1,:) + 1.2665468516202070398124490679947*x(1,:)^2 + 0.058711157666221952533547323582752*x(1,:)^3 - 0.28781069452828345056616399233462*x(1,:)^4 - 0.29722062523126249744542803910008;

% dxdt(1,:) = x(3,:)^2+x(2,:);
% dxdt(2,:) = -0.07225017783189037324910310012961*x(1,:)^3-1.0000000000000000176603214988472*x(1,:)^2+0.28590872180264067144624580881403*x(1,:)+x(3,:)+0.00000000000000011226827621331125672584339169935;
% dxdt(3,:) = -x(1,:)^2*x(3,:)+4.2408902470083598146288750285748*x(3,:)^4-0.0000000000000029785537508371308782725154515509*x(3,:)^3-8.12955283323774979820086628024*x(3,:)^2+0.0000000000000016462861212869430297641574498419*x(3,:)+1.9543387965142491324854745471384;

%% 3-3
% dxdt(1,:) = x(3,:)^2+x(2,:);
% dxdt(2,:) = -0.092329208875508065282836603863096*x(1,:)^3-1.0000000000000000034874012454557*x(1,:)^2+0.84921676986354214614986328039474*x(1,:)+x(3,:)-0.000000000000000031971192824431822723110163777761;
% dxdt(3,:) = -x(1,:)^2*x(3,:)-4.9857772792774355252731766086072*x(3,:)^3-0.000000000000000062773222418202686345491548390287*x(3,:)^2+5.0953006191812528768991796823684*x(3,:)-0.000000000000000063942385648863645446220327555522;

%% sos-model 3-4
dxdt(1,:) = 0.0000000000016975534264690627264748450605129*x(1,:)^4+0.02698227887897345315074651942403*x(1,:)^3-1.0000000000089672340586282760011*x(1,:)^2-0.23073553621007529083423529906819*x(1,:)+1.4910341041517696123497198563345;
dxdt(2,:) = -x(2,:)*x(1,:)^3-x(2,:);
dxdt(3,:) = -x(1,:)^2*x(3,:)-3.1534768207698271602623663056875*x(3,:)^4-3.8452229138132887342749199888203*x(3,:)^3+1.7197670368998982937114305968862*x(3,:)^2+1.0541896624154754036339909362141*x(3,:)-1.1253259617261148761713229760062;

end