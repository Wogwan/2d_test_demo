clear;
tic
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
f = [x2-x1;
    0.46382156330341729589940137480476*x1^4-0.21219821363526560116570580989732*x1^3+0.35000677410250802565647474138031*x1^2+x1*x2-0.27258976376810471290805063896793*x1+0.15407267311901634565529661813343*x2^4+0.94268273772394006737584959410015*x2^3+4.1277331008261466394060335005634*x2^2-1.7098389065959235244562819389103*x2+0.011510115809394316777058975276304
    ];
gg = [1;1];
input = [gg(1)*u1;gg(2)*u2];
%%
% f = [x2; -x1-x2*(1-x1^2)];
% gg = [1;1];
% input = [gg(1)*u1;gg(2)*u2];
% V = x1^2+x1*x2+x2^2;
V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
%%
% V = x1^2+x1*x2+x2^2;
% V = 1*x1^4+1*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
%%
C0 = 0.1;
cc = 1;
k = 2;
k_l = 4;
dom = 10;
domain = [-dom dom -dom dom];
%%
C1 = (x1+4)^2+(x2-6)^2-4;
C2 = (x1+3)^2+(x2+4)^2-4;
C3 = (x1-6)^2+(x2-0)^2-5;
C = [C1;C2;C3];
% C = [C1;C2];
kk = 1;
%%
figure_id = 11;
figure(figure_id);clf;hold on;
[~,~]=pcontour(C(1),0,domain,'r');            % Plot the original Lyapunov sublevel set
[~,~]=pcontour(C(2),0,domain,'r');            % Plot the original Lyapunov sublevel set
[~,~]=pcontour(V,C0,domain,'g');              % Plot the original Lyapunov sublevel set
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'r');            % Plot the original Lyapunov sublevel set
end

solU = [];
v_c = [];
iter = 1;
%%
while abs(double(cc)-double(C0)) >= 1e-8
    iter = iter + 1;
    if iter ~= 1
        C0 = cc;
    end
%     [solu1,solu2,solL,kk]= sos_function_v(f,gg,k,k_l,V,C0);
    [solu1,solu2,solL,kk] = sos_function_v(f,gg,k,k_l,V,C0);
    if kk == 0
        break
    end
%     [cc,kk,solu1,solu2] = sos_function_v2(f,gg,k,k_l,V,C,dom,solL);
    [cc,kk,solu1,solu2] = sos_function_v2(f,gg,k,k_l,V,C,dom,solL,figure_id)
    v_c = [v_c; double(cc)]
    solU = [solU;[solu1,solu2]];
    if kk == 0
        figure(figure_id);hold on;
        [~,~]=pcontour(V,v_c(end),domain,'b');  
        break
    end
end
toc