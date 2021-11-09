clear;close all;
diary off
diary Test_check_system_1.out;
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%% Initial System information
C1 = (x1+8)^2+(x2-0)^2-4;
C2 = (x1-8)^2+(x2+0)^2-4;
C3 = (x1-0)^2+(x2-8)^2-4;
C4 = (x1-0)^2+(x2+8)^2-4;
% sym.us = [C1;C2;C3;C4];
sym.us = [C1;C2;C3];
sym.f = [x2-x1;
    0.46382156330341729589940137480476*x1^4-0.21219821363526560116570580989732*x1^3+0.35000677410250802565647474138031*x1^2+x1*x2-0.27258976376810471290805063896793*x1+0.15407267311901634565529661813343*x2^4+0.94268273772394006737584959410015*x2^3+4.1277331008261466394060335005634*x2^2-1.7098389065959235244562819389103*x2+0.011510115809394316777058975276304
    ];
% sym.f = [x1^2*x2^2-2*x1-x2;
%     0.050980740770243615500589839939494*x1^4-1.2396279563068186568841611006064*x1^3-x1^2*x2+5.6603625282381084815597205306403*x1^2+4.2997830226381255069867393103777*x1+0.007305173261451977997915641083182*x2^4+0.010718334556742890872893525511245*x2^3-0.12879414181748349843559253713465*x2^2+0.9648658043110491799865258144564*x2-0.024152264962771043121936287434437];
sym.gg = [1; 1]; sym.sys = [sym.f(1)+sym.gg(1)*u1; sym.f(2)+sym.gg(2)*u2];
sym.V0 = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
% b = p2s([x1^4, x1^2*x2^2, x1^2, x1*x2, x2^4, x2^2]);
% llp1_vector = [-5 -4 -3 -2 -1 1 2 3 4 5 6];
% llp2_vector = [-5 -4 -3 -2 -1 1 2 3 4 5 6];
% llp3_vector = [-5 -4 -3 -2 -1 1 2 3 4 5 6];
% llp4_vector = [-5 -4 -3 -2 -1 1 2 3 4 5 6];
% llp5_vector = [-5 -4 -3 -2 -1 1 2 3 4 5 6];
% llp6_vector = [-5 -4 -3 -2 -1 1 2 3 4 5 6];


iter_id = 1;

% for i1 = 1:length(llp1_vector(1,:))
%     lp1 = llp1_vector(1,i1);
%     for i2 = 1:length(llp2_vector(1,:))
%         lp2 = llp2_vector(1,i2);
%         for i3 = 1:length(llp3_vector(1,:))
%             lp3 = llp3_vector(1,i3);
%             for i4 = 1:length(llp4_vector(1,:))
%                 lp4 = llp4_vector(1,i4);
%                 for i5 = 1:length(llp5_vector(1,:))
%                     lp5 = llp5_vector(1,i5);
%                     for i6 = 1:length(llp5_vector(1,:))
%                         lp6 = llp6_vector(1,i6);
%                         iter_id = iter_id + 1
%                         coe = [lp1;lp2;lp3;lp4;lp5;lp6];
%                         mid = b*coe;
%                         sym.V0 = s2p(mid);

%% lp1.v_k_u = 2; % lp1.v_k_l = 4;
lp1_vector = [2 4;4 4];

%% lp2.u = 4; % lp2.h = 4; % lp2.us = 8; % lp2.au = 8; % lp2.gamma = 0;
lp2_vector = [2 4 4 4 0;2 4 4 4 1;4 4 4 4 0];

%% lp3.V_us = 6; % lp3.V_au = 8; % lp3.V_degree = 4; % lp3.gamma = 0;
lp3_vector = [6 8 4 0; 4 4 4 0];

%% lp4.u = 4;% lp4.au = 8;
lp4_vector = [4 4;4 6;4 8;6 4;6 6;6 8;8 8];

%% lp5.u = 4; % lp5.h = 4; % lp5.us = 8; % lp5.au = 8; % lp5.gamma = 0;
lp5_vector = [4 4 8 8 0;
    4 4 8 8 0];

%%
iter = 1;
plot.figure_id = 1;

%%
for i1 = 1:length(lp1_vector(:,1))
    lp1 = struct('v_k_u',lp1_vector(i1,1),'v_k_l',lp1_vector(i1,2));
    for i2 = 1:length(lp2_vector(:,1))
        lp2 = struct('u',lp2_vector(i2,1),'h',lp2_vector(i2,2),'us',lp2_vector(i2,3),'au',lp2_vector(i2,4),'gamma',lp2_vector(i2,5));
        for i3 = 1:length(lp3_vector(:,1))
            lp3 = struct('V_us',lp3_vector(i3,1),'V_au',lp3_vector(i3,2),'V_degree',lp3_vector(i3,3),'gamma',lp3_vector(i3,4));
            for i4 = 1:length(lp4_vector(:,1))
                lp4 = struct('u',lp4_vector(i4,1),'au',lp4_vector(i4,2));
                for i5 = 1:length(lp5_vector(:,1))
                    iter = iter + 1
                    lp5 = struct('u',lp5_vector(i5,1),'h',lp5_vector(i5,2),'us',lp5_vector(i5,3),'au',lp5_vector(i5,4),'gamma',lp5_vector(i5,5));
                    kkk = auto_test_demo_2d(sym,lp1,lp2,lp3,lp4,lp5,plot)
                    plot.figure_id = plot.figure_id + kkk;
                end
            end
        end
        
    end
end
%                     end
%                 end
%             end
%         end
%     end
%
% end
diary off