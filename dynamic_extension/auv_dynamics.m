function dxdt = auv_dynamics(t,x)
% Differential equations for the AUV
    c=0.02567;
    a1=-0.28;a2=-0.993;a3=-0.1729;
    b1=0.0538;b2=-0.0122;b3=-0.0349;
    
    % Predefined initial positions, final positions and velocity
    global p1 p2 q1 q2 v

    % Defining desired trajectory (a straight line in this case)
    x1d=p1+v*(q1-p1)*t;
    x2d=p2+v*(q2-p2)*t;
    x1d_dot=v*(q1-p1);
    x2d_dot=v*(q2-p2);
    x1d_ddot=0;
    x2d_ddot=0;
    x1d_dddot=0;
    x2d_dddot=0;

    % Testing for arbitrary desired trajectory
%     x1d=t/10;
%     x2d=sin(t/10);
%     x1d_dot=0.1;
%     x2d_dot=0.1*cos(t/10);
%     x1d_ddot=0;
%     x2d_ddot=-0.01*sin(t/10);
%     x1d_dddot=0;
%     x2d_dddot=-0.0001*cos(t/10);

    % Defining states of the transformed system 
    x1=x(1);
    x1_dot=x(4)*cos(x(3))-x(5)*sin(x(3));
    x1_ddot=(a1*x(4)+b1*x(5)*x(6)+x(7))*cos(x(3))-c*x(4)*x(6)*sin(x(3))...
        -(a2*x(5)+b2*x(4)*x(6))*sin(x(3))-c*x(5)*x(6)*cos(x(3));
    x2=x(2);
    x2_dot=x(4)*sin(x(3))+x(5)*cos(x(3));
    x2_ddot=(a1*x(4)+b1*x(5)*x(6)+x(7))*sin(x(3))+c*x(4)*x(6)*cos(x(3))...
        +(a2*x(5)+b2*x(4)*x(6))*cos(x(3))-c*x(5)*x(6)*sin(x(3));

    eval1=0.2;eval2=0.3;eval3=0.4;
    eval4=0.2;eval5=0.3;eval6=0.4;
    
    % Defining desired control input for asymptotic tracking
    % of error in x-y positions
    ux=x1d_dddot+(eval1+eval2+eval3)*(x1d_ddot-x1_ddot)...
        +(eval1*eval2+eval1*eval3+eval2*eval3)*(x1d_dot-x1_dot)+(eval1*eval2*eval3)*(x1d-x1);
    uy=x2d_dddot+(eval4+eval5+eval6)*(x2d_ddot-x2_ddot)...
        +(eval4*eval5+eval4*eval6+eval5*eval6)*(x2d_dot-x2_dot)+(eval1*eval2*eval3)*(x2d-x2);
    
    f1=a1^2*x(4)+((a1+a2+a3)*b1-2*c*a2-c*a3)*x(5)*x(6)...
        +(b1*b2-2*c*b2-c^2)*x(4)*x(6)^2+(b1-c)*b3*x(4)*x(5)^2;
    f2=a2^2*x(5)+(a1+a2+a3)*b2*x(4)*x(6)+(b1*b2+2*c*b1-c^2)*x(5)*x(6)^2 ...
        +(b2*b3+c*b3)*x(4)^2*x(5)+(2*c*a1+c*a3)*x(4)*x(6);

    tau2=(-(b2+2*c)*x(6)*x(7)-f2-ux*sin(x(3))+uy*cos(x(3)))/((b2+c)*x(4));
    
    % Dynamics
    dxdt=[x(4)*cos(x(3))-x(5)*sin(x(3));
        x(4)*sin(x(3))+x(5)*cos(x(3));
        c*x(6);
        a1*x(4)+b1*x(5)*x(6)+x(7);
        a2*x(5)+b2*x(4)*x(6);
        a3*x(6)+b3*x(4)*x(5)+tau2;
        -a1*x(7)-(((b1-c)*x(5))/((b2+c)*x(4)))*(-(b2+2*c)*x(6)*x(7)-f2-ux*sin(x(3))+uy*cos(x(3)))-f1+ux*cos(x(3))+uy*sin(x(3))];
    
end