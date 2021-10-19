function solu = sos_function_3(k,h,gamma,mm,V,C,V0)
    pvar x1 x2;
    x = [x1;x2];
    f = [sys2pvar(dotx2(1));
        sys2pvar(dotx2(2))];                                      % Convert the doxt3 from syms form to pvar form for SOS analysic
    V = 1*(2*x1^4+0.5*x2^4-x1^2*x2^2)+1*(1*x1^2+1*x2^2-0.5*x1*x2);

    episC = 1e-4;                                                 % The condition to stop searching sub-level set C         
    epsi = 1e-5;                                                  % The condition to stop searching barrier certificate
    C0 = 1e-2;                                                    % The initial value for searching sub-level set C
    C = 0.1;                                                      % The initial value for searching sub-level set C
    k = 0;                                                        % The initial value for recording steps for sub-level set C
    %     domain = [-3 3 -3 3];
%     domain = [-1.5 1.5 -1.5 1.5];                                 % The condition to stop searching barrier certificate
    C_L_factor = 4;

    [~,~]=pcontour(V,double(C0),domain,'c'); hold on;             % Plot the original Lyapunov sublevel set
    Vdot = jacobian(V,x)*f;
    [~,~]=pcontour(Vdot,0,domain,'k'); hold on;
    
    while double(C)-double(C0)>=episC
        k = k + 1;
        fprintf('K =   %d\n ',k);
        if k==0
            L = search_LC_pvar_main(f,V,C0,C_L_factor);
        else
            C0 = C;
            [L, k_mark] = search_LC_pvar_main(f,V,double(C),C_L_factor);
        end
        if k_mark == 1 && k == 1
            break
        end
        C = search_C_pvar_main(f,V,L);
        refreshdata; drawnow; 
    end
    
end