pvar x1 x2 u htol epsi;
x = [x1;x2];

f = [0.1*x1^2+1*x2; 0.1*x1*x2-0.2*x1+(1+x1^2)*u];
V = 1*x1^2+1*x2^2+1*x1*x2;
C0 = 3;
cc = 3.2;
dom = 4;
domain = [-dom dom -dom dom];
C1 = (x1-3)^2+(x2-1)^2-1;
C2 = (x1+3)^2+(x2+4)^2-1;
C3 = (x1+4)^2+(x2-5)^2-1;
C4 = (x1-2)^2+(x2-6)^2-1;
C5 = -x2+2;
C6 = (x1+0.5)^2+(x2+3.25)^2-1;
C7 = -x2 - 4;
C8 = -x1 + 1;
C = [C1;C2;C3;C4;C5;C6;C7;C8];

figure(11);clf;
[~,~]=pcontour(C(1),0,domain,'r'); hold on;             % Plot the original Lyapunov sublevel set

k = 2;
iter = 1;
while double(cc)-double(C0) >= 1e-3
    iter = iter + 1;
    if iter ~= 1
        C0 = cc;
    end
    [solu,solL,kk]=sos_function_v(f,k,V,C0,dom);
    [cc,kk]=sos_function_v2(f,k,V,C,dom,solL);
end