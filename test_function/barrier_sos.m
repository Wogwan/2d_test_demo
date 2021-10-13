pvar x1 x2 gam;
x = [x1;x2];

deg = 4;
k = deg;
dom = 4;
domain = [-dom dom -dom dom];

% Create vector field: dx/dt = f(x)/(x3^2+1)
f = [x2; 
    (1-x1^2)*x2-x1];

gam_0 = 0.4;
c2 = gam_0 -(x1^2+x2^2); % Safe region qualified polynomial X_0
c3 = x1^2+x2^2-4;      % Unsafe region : X_u
c4 = x1^2+x2^2;        % For all [x x] : X
c5 = c4-c3;            % X \ X_u

[~,~]=pcontour(c2,0,domain,'c'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(c3,0,domain,'c'); hold on;

episC = 1e-8;   
gam = 0.5;
k = 0;

while abs(double(gam)-double(gam_0))>=episC
    fprintf('K =   %d\n ',k+1);
    if k==0
        [Bs, S] = barrier_test(deg,f,double(gam_0),[c2,c3,c4,c5]);
    else
        gam_0 = gam;
        c2 = gam_0-c4;
        [Bs, S] = barrier_test(deg,f,double(gam_0),[c2,c3,c4,c5]);
    end
    [B, gam] = barrier_gam(deg,f,B,S,[c2,c3,c4,c5]);
end

%%
plotVF_validate;hold on;
[~,~]=pcontour(B,double(gam),domain,'b'); hold on;             % Plot the original Lyapunov sublevel set
figure(2);hold on;
xlim([-dom dom]); ylim([-dom dom]); hold on;   
[~,~]=pcontour(c2,0,domain,'m'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(c3,0,domain,'m'); hold on;             % Plot the original Lyapunov sublevel set
[~,~]=pcontour(Bs,0,domain,'g'); hold on;             % Plot the original Lyapunov sublevel set
Bdot = jacobian(Bs,x)*f;
[~,~]=pcontour(Bdot,0,domain,'b'); hold on;
axis(domain);
