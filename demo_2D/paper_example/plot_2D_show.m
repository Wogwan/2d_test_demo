% clear;
pvar x1 x2;
%%
% load BC_test_2D_after_meeting;
%%
dom = 10; domain = [-dom dom -dom dom];
C1 = (x1+4)^2+(x2-6)^2-4;
C2 = (x1+3)^2+(x2+4)^2-4;
C3 = (x1-6)^2+(x2-0)^2-5;
%%
C = [C1;C2;C3];
k = ['r','g','b','m','c','k','y'];
for iter = 1:length(A)
    %     B = A(iter,3);
    B = Barrier(iter)
    if mod(iter,30) == 0
        pause
    end
    figure(1);clf;hold on;
    [~,~]=pcontour(C(1),0,domain,'r');
    [~,~]=pcontour(C(2),0,domain,'r');
    [~,~]=pcontour(C(3),0,domain,'r');
    if mod(iter,7) == 0
        [~,h31]=pcontour(B,0,domain,k(7)); hold on;
        h31.LineStyle = '-'; h31.LineWidth = 1.4;
    else
        [~,h31]=pcontour(B,0,domain,k(mod(iter,7))); hold on;
        h31.LineStyle = '-'; h31.LineWidth = 1.4;
    end
    xlim([-dom dom]); ylim([-dom dom]);
end

