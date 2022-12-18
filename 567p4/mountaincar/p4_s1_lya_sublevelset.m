function sys = p4_s1_lya_sublevelset(sys,Lya_index)
%%
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
f = s2p(sys.f);
gg = sys.gg;
V = sys.V(Lya_index);
C = sys.us_region;
dom = sys.dom;
domain = sys.domain;
%%
C0 = 0.1; cc = 1;
k = 4; k_l = 4;
kk = 1; iter = 1;
figure_id = 11;
solU = []; v_c = [];
%%
while double(cc)-double(C0) >= 1e-6
    iter = iter + 1;
    if iter ~= 1
        C0 = cc;
    end
    [solu1,solu2,solL,kk] = p4_sos_function_v1(f,gg,k,k_l,V,C0);
    if kk == 0
        break
    end
    [cc,kk,solu1,solu2] = p4_sos_function_v2(f,gg,k,k_l,V,C,dom,solL,figure_id);
    v_c = [v_c; double(cc)];
    solU = [solU;[solu1,solu2]];
end
if sys.plott == 1
    figure(figure_id);clf;hold on;
    [~,~]=pcontour(V,v_c(end),domain,'b');
end
%%
% fprintf('Sublevel set of V0(x) is %.14f. \n',v_c(end));
sys.cv = double(v_c(end));
if v_c(end)==0
   sys.exist_bar = 0; 
end
sys.c_U = solU(end,2);
end