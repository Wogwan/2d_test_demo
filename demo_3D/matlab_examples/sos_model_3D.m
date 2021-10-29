function [dx,f,g] = sos_model_3D(k)
syms x1 x2 x3;
x = [x1 x2 x3];
switch k
    case '3d_1'
        f = [x(2)+x(3)^2;
            x(3)-x(1)^2;
            -x(1)-2*x(2)-x(3)+x(2)^3];
        g = [0;
            -x(1)*sin(x(1));
            -(exp(2*x(3))*cos(x(3))^2)^(1/4)];
        dx = f+g;
    case '3d_2'
        f = [x(2)+x(3)^2;
            x(3)-x(1)^2;
            -x(1)^2*x(3)];
        g = [0;
            1*sin(2*x(1));
            2*cos(x(3))];
        dx = f+g;
    case '3d_3'
        f = [x(2)+x(3)^2;
            x(3)-x(1)^2;
            -x(1)^2*x(3)];
        g = [0;
            1*sin(1*x(1));
            2*sin(x(3))];
        dx = f+g;
    case '3d_4'
        f = [-x(1)*x(3)^2;
            -x(1)*x(3);
            -x(1)^2*x(3)];
        g = [0;
            1-x(1)*cos(x(2))^2;
            0];
        dx = f+g;
    otherwise
        f = [ -x(1)-2*x(2)+x(1)*(4*x(1)^2-2*x(1)*x(2)+x(2)^2);
            2*x(1)+x(2)*(4*x(1)^2-2*x(1)*x(2)+x(2)^2)];
        g = [0;
            sqrt(sqrt((exp(x(:,1))*cos(x(:,1)))^2))];
        dx = f+g;
end
end