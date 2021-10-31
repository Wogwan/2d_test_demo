clear;
tic
pvar x1 x2 u1 u2;
x = [x1;x2];
%%
f = [x2-x1
    0.15432387994108778662817450645485*x1^4+0.23541948636981062330190921043309*x1^3+0.41058519554981867980300888929277*x1^2+x1*x2-0.46368439821849470143059572061854*x1+0.018934809918422834673634724822477*x2^4+0.053017123265488262651157214122577*x2^3+0.1081959091461522914912052328873*x2^2-0.9986920403516159766305754219573*x2-0.0018682553605258930828902919074608
    ];
gg = [1;1];
input = [gg(1)*u1;gg(2)*u2];
%%
C0 = 100;
cc = 11;
k = 2;
k_l = 2;
L_au = 2;
dom = 10;
domain = [-dom dom -dom dom];
%%
C1 = (x1+2)^2+(x2-3)^2-2;
C2 = (x1+3)^2+(x2+2)^2-1;
C3 = (x1-2)^2+(x2-2)^2-1;
%%
C4 = (x1-2)^2+(x2-6)^2-1;
C5 = -x2+2;
C6 = (x1+0.5)^2+(x2+3.25)^2-1;
C7 = -x2 - 4;
C8 = -x1 + 1;
C = [C1;C2;C3;C4;C5;C6;C7;C8];
kk = 1;
%%
B = -20920.02735971049*x1^2-12944.00205666696*x1*x2-21004.89487163851*x2^2+725.8021517955463*x1-34941.78371754049*x2+132960.4633676535;
u = -1.423136761279898*x1^2+0.1389802653220975*x1*x2+0.1251182246445742*x2^2-1.103317535267957*x1-1.214738774413024*x2-0.02619321979348357;
%%
figure(11);clf;hold on;
[~,~]=pcontour(C(1),0,domain,'r');            % Plot the original Lyapunov sublevel set
[~,~]=pcontour(C(2),0,domain,'r');            % Plot the original Lyapunov sublevel set
[~,~]=pcontour(C(3),0,domain,'r');            % Plot the original Lyapunov sublevel set
[~,~]=pcontour(B,0,domain,'b');            % Plot the original Lyapunov sublevel set
hold on;
%%
[V,kk]=sos_optimal_v(f,gg,k_l,B,u,C);
% [~,~]=pcontour(V,C0,domain,'g');              % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(V,0,domain,'g');              % Plot the original Lyapunov sublevel set

%%
solU = [];
v_c = [];
iter = 1;

%%
% while abs(double(cc)-double(C0)) >= 1e-7
%     iter = iter + 1;
%     if iter ~= 1
%         C0 = cc;
%     end
%     [solu,solL,kk]= sos_function_v(f,gg,k,k_l,V,C0,dom);
%     [cc,kk,solu] = sos_function_v2(f,gg,k,k_l,V,C,dom,solL);
%     v_c = [v_c; double(cc)]
%     solU = [solU;solu];
%     if kk == 0
%         figure(11);hold on;
%         [~,~]=pcontour(V,v_c(end),domain,'g');  
%         break
%     end
% end

%%
Control = [];
% solh = V;
[a1,b1] = coeffs(p2s(V));
C0 = vpa(a1(end));
b1(end)=[];
a1(end)=[];
V = s2p(a1*b1');
[~,~]=pcontour(V,0,domain,'g');
% solh = s2p(a1*b1');
mm = 0;
gamma = 0;
solh = s2p(C0) + V;
[~,~]=pcontour(V,C0,domain,'R');              % Plot the original Lyapunov sublevel set
% [~,~]=pcontour(solh,0,domain,'b');              % Plot the original Lyapunov sublevel set


%%
for i = 1:40
    fprintf('i=%6.0f\n',i);
%     [SOLu,SOL1,SOL2,kk] = sos_optimal_v2(f,gg,k,L_au,B,C0);
    [SOLu,SOL1,SOL2,kk] = sos_function_1(f,k,solh,V,mm,gamma,gg,L_au);
    Control = [Control; SOLu];
    if kk == 0
        break
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [solh,trace_Q,kk]=sos_function_2(f,k,SOLu,SOL1,SOL2,gamma,mm,V,C,dom,gg,L_unsafe_factor);
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
    if kk == 0
        break
    end
    %%%%%%%
end

toc