function [h, kk] = sos_optimal_V1_2D(f,gg,B,u1,u2,l_au,l_us,h_degree,C,gamma,V)

pvar x1 x2 Vtol2 cc1 cc2 cc3;
x = [x1;x2];
%%
% [V,vc] = polydecvar('v_w',[x1^4; x2^4; x1^2*x2^2;x1^2;x2^2;x1*x2]);
[h,hc] = polydecvar('h_w',monomials(x,0:h_degree));
%%
% [L1,L1_Q] = polydecvar('L1_w',monomials(x,0:l_au));
% [L2,L2_Q] = polydecvar('L2_w',monomials(x,0:l_au)); 
% [L3,L3_Q] = polydecvar('L3_w',monomials(x,0:l_us)); 
% [L4,L4_Q] = polydecvar('L4_w',monomials(x,0:l_us)); 
% [L5,L5_Q] = polydecvar('L5_w',monomials(x,0:l_us)); 
% % [L6,L6_Q] = sosdecvar('L6_w',monomials(x,0:l_us/2)); 
%%
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:l_au/2));
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:l_au/2)); 
% [L2,L2_Q] = polydecvar('L2_w',monomials(x,0:l_au/2)); 
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:l_us/2)); 
[L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:l_us/2)); 
[L5,L5_Q] = sosdecvar('L5_w',monomials(x,0:l_us/2)); 
[L6,L6_Q] = sosdecvar('L6_w',monomials(x,0:l_au/2)); 
%%
hdot = jacobian(h, x1)*(f(1)+gg(1)*u1)+ jacobian(h, x2)*(f(2)+gg(2)*u2);
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
%% Constraint
pcr_11 = L1 >= 0;
% pcr_12 = L2 >= 0;
pcr_13 = L3 >= 0;
pcr_14 = L4 >= 0;
pcr_15 = L5 >= 0;
pcr_16 = L6 >= 0;
%%
pconstr_1 = h-L1*B-Vtol2 >= 0;
% pconstr_1 = -V+L1*B+Vtol2 <= 0;
pconstr_2 = -hdot-L2*B-gamma*B+Vtol2 <= 0;
% pconstr_2 = Vdot-hdot-L2*B-gamma*B+Vtol2 <= 0;
pconstr_3 = Vtol2 >= 0;
% pconstr_4 = -hdot-L6*Vdot >= 0;
pconstr_4 = -hdot-L6*Vdot >= 0;
%%
pcr_21 = h-C(1)*L3 <= 0;
pcr_22 = h-C(2)*L4 <= 0;
if length(C) == 2
    pconstr = [pcr_11;pcr_12;pcr_13;pcr_14;pconstr_1;pconstr_2;pconstr_3;pcr_21;pcr_22];
elseif length(C) == 3
    pcr_23 = h-C(3)*L5 <= 0;
%     pconstr = [pcr_11;pcr_12;pcr_13;pcr_14;pcr_15;pconstr_1;pconstr_2;pconstr_3;pcr_21;pcr_22;pcr_23];
%     pconstr = [pcr_11;pcr_12;pcr_13;pcr_14;pcr_15;pcr_16;pconstr_1;pconstr_2;pconstr_3;pconstr_4;pcr_21;pcr_22;pcr_23];
%     pconstr = [pcr_11;pcr_13;pcr_14;pcr_15;pconstr_1;pconstr_2;pconstr_3;pcr_21;pcr_22;pcr_23];

pconstr = [pcr_11;pcr_13;pcr_14;pcr_15;pcr_16;pconstr_1;pconstr_2;pconstr_3;pconstr_4;pcr_21;pcr_22;pcr_23];

else
    fprintf('Constraints vector does not match.======\n');
end  

%% Inverse constraint
% pcr_21 = V+C(1)*L3 >= 0;
% pcr_22 = V+C(2)*L4 >= 0;
% if length(C) == 2
%     pconstr = [pcr_11;pcr_12;pcr_13;pcr_14;pconstr_1;pconstr_2;pconstr_3;pcr_21;pcr_22];
% elseif length(C) == 3
%     pcr_23 = V+C(3)*L5 >= 0;
%     pconstr = [pcr_11;pcr_12;pcr_13;pcr_14;pcr_15;pconstr_1;pconstr_2;pconstr_3;pcr_21;pcr_22;pcr_23];
% else
%     fprintf('Constraints vector does not match.======\n');
% end  

%%
obj = Vtol2;
% obj = -Vtol2;
%%
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,obj,opts);
if info.feas
    kk = 1;
    h = subs(h,dopt)
else
    kk = 0;
    h  = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end

end