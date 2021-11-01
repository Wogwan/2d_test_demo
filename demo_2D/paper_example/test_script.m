pvar x1 x2;
x = [x1;x2];
dom = 1000; domain = [-dom dom -dom dom];
% V = -0.0002688249895716166*x1^2+4.782656596144253e-05*x1*x2+0.001035911584956916*x2^2+6.53143778574437e-05*x1+9.265590495178159e-05*x2+722.2262195419746;
figure(15);clf;hold on;
cc = [1,5,10];
for i = 1:10
    hold on;
    [~,~]=pcontour(V,cc(1)/(i^1),domain,'r');
    [~,~]=pcontour(V,cc(2)/(i^1),domain,'g');
    [~,~]=pcontour(V,cc(3)/(i^1),domain,'b');
end