clear;
pvar x1 x2;
%%
% load Opt_Ly_2_BROA;
% load BC_test_2D;
%%
% load Opt_Ly_2_BROA2;
% load BC_test_2D2;
load BC_test_2D_after_night_train
%%
dom = 10; domain = [-dom dom -dom dom];
iter = 1;
% C1 = (x1+5)^2+(x2-8)^2-6;
% C2 = (x1+7)^2+(x2+7)^2-6;
% C3 = (x1-6)^2+(x2-0)^2-6;
%%
C1 = (x1-4)^2+(x2+4)^2-3;
C2 = (x1+4)^2+(x2+6)^2-6;
C3 = (x1-3)^2+(x2-5)^2-2;
%%
C = [C1;C2;C3];
k = ['r','g','b','m','c','k','y'];
for iter = 1:length(A)
%     V0 = V_Sel(iter);
%     C0 = C0_Sel(iter);
    B = A(iter,3);
    B_ad = A(iter,4);
    if mod(iter,30) == 0
        pause
    end
    figure(1);clf;hold on;
    [~,~]=pcontour(C(1),0,domain,'r');
    [~,~]=pcontour(C(2),0,domain,'r');
    [~,~]=pcontour(C(3),0,domain,'r');
    if mod(iter,7) == 0
        [~,~]=pcontour(B,0,domain,k(7)); hold on;             % Plot the original Lyapunov sublevel set
        [~,h31]=pcontour(B_ad,0,domain,'m'); hold on;             % Plot the original Lyapunov sublevel set
        h31.LineStyle = '-'; h31.LineWidth = 1.4;
    else
        [~,~]=pcontour(B,0,domain,k(mod(iter,7))); hold on;             % Plot the original Lyapunov sublevel set
        [~,h31]=pcontour(B_ad,0,domain,'m'); hold on;             % Plot the original Lyapunov sublevel set
        h31.LineStyle = '-'; h31.LineWidth = 1.4;
    end
%     [~,h32]=pcontour(V0,double(C0),domain,'k'); hold on;             % Plot the original Lyapunov sublevel set
%     h32.LineStyle = '-'; h32.LineWidth = 2;
    xlim([-dom dom]); ylim([-dom dom]);
    %     pause
end

