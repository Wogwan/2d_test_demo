% pvar x1 x2;
% x = [x1;x2];
dom = 20; domain = [-dom dom -dom dom];
xlim([-dom dom]);ylim([-dom dom]);
% V = -0.0002688249895716166*x1^2+4.782656596144253e-05*x1*x2+0.001035911584956916*x2^2+6.53143778574437e-05*x1+9.265590495178159e-05*x2+722.2262195419746;
figure(15);clf;hold on;
cc = [1,2,3];
for i = 0:16
    k = ['r','g','b','m','c','k','y'];
    if mod(i,7) == 0
        [~,~]=pcontour(V,cc(1)*i,domain,k(7)); hold on;
        [~,~]=pcontour(V,cc(2)*i,domain,k(7)); hold on;
        [~,~]=pcontour(V,cc(3)*i,domain,k(7)); hold on;             
    else
        [~,~]=pcontour(V,cc(1)*i,domain,k(mod(i,7))); hold on;      
        [~,~]=pcontour(V,cc(2)*i,domain,k(mod(i,7))); hold on;      
        [~,~]=pcontour(V,cc(3)*i,domain,k(mod(i,7))); hold on;      
    end
end
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
figure(16);clf;dom = 100;
for i = 0:16
    k = ['r','g','b','m','c','k','y'];
    if mod(i,7) == 0
        [~,~]=pcontour(Vdot,0,domain,k(7)); hold on;             % Plot the original Lyapunov sublevel set
    else
        [~,~]=pcontour(Vdot,0,domain,k(mod(i,7))); hold on;             % Plot the original Lyapunov sublevel set
    end
end