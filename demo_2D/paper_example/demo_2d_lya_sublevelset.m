clear;
tic
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
f = [x2-x1
    -0.23606416828637188828796009029584*x1^4+0.24650838581310801777701368775053*x1^3+x1^2*x2-0.7173309565565305634393666878168*x1^2+0.26453823849708858750862106035129*x1+0.03729676680157056195552556232542*x2^4+0.01893812978180069162004173222158*x2^3+0.13680422363948308017711497086566*x2^2+0.043563863537901807709840085180986*x2+0.23883336887991710173473336453753
    ];
gg = [1;1];
input = [gg(1)*u1;gg(2)*u2];
% V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; 
V = 1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1;
%%
C0 = 0.1;
cc = 1;
k = 4;
k_l = 4;
dom = 15;
domain = [-dom dom -dom dom];
%%
C1 = (x1+4)^2+(x2-6)^2-4;
C2 = (x1+0)^2+(x2+5)^2-4;
C3 = (x1-5)^2+(x2-0)^2-4;
%%
C = [C1;C2;C3];
kk = 1;
%%
figure_id = 11;
figure(figure_id);clf;hold on;
[~,~]=pcontour(C(1),0,domain,'r');          
[~,~]=pcontour(C(2),0,domain,'r');          
[~,~]=pcontour(V,C0,domain,'g');            
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'r');            
end
%%
solU = [];
v_c = [];
iter = 1;
%%
while double(cc)-double(C0) >= 1e-6
    iter = iter + 1;
    if iter ~= 1
        C0 = cc;
    end
    [solu1,solu2,solL,kk] = sos_function_v(f,gg,k,k_l,V,C0);
    if kk == 0
        break
    end
    [cc,kk,solu1,solu2] = sos_function_v2(f,gg,k,k_l,V,C,dom,solL,figure_id);
    v_c = [v_c; double(cc)]
    solU = [solU;[solu1,solu2]];
    if kk == 0
        figure(figure_id);hold on;
        [~,~]=pcontour(V,v_c(end),domain,'b');
        break
    end
end
figure(figure_id);hold on;
[~,~]=pcontour(V,v_c(end),domain,'b');
%%
%%
fprintf('Sublevel set of V0(x) is %.14f. \n',v_c(end));
toc