function dxdt = auv_dynamics(t,x)
% Differential equations for the AUV
    c=0.02567;
    a1=0.28;a2=0.993;a3=0.1729;
    b1=0.0538;b2=0.0122;b3=0.0349;
    y=coordinate_transform(x);
    
    % For bringing to origin
%     v1=-2*y(1)-3*y(2);
%     v2=-2*y(3)-3*y(4); 

    % For tracking    

    % First reference trajectory
    
%     p1=1;p2=1;
%     q1=2;q2=3;
%     v=1;
%     x1d=p1+v*(q1-p1)*t;
%     x2d=p2+v*(q2-p2)*t;
%     x3d=atan2(x2d,x1d);
%     x1d_dot=v*(q1-p1);
%     x2d_dot=v*(q2-p2);
%     x1d_ddot=0;
%     x3d_dot=(v*(p1*q2-q1*p2))/(x1d^2+x2d^2);
%     x3d_ddot=((v*(-p1*q2+q1*p2))/(x1d^2+x2d^2)^2)*(2*x1d*x1d_dot+2*x2d*x2d_dot);
%     v1=x1d_ddot-3*(y(2)-x1d_dot)-2*(y(1)-x1d);
%     v2=x3d_ddot-3*(y(4)-x3d_dot)-2*(y(3)-x3d);
    
    % Second reference trajectory
    % Another edit is possible if we consider only the slope and not the
    % final points, maybe that way there will be no rush for the AUV to
    % track the actual trajectory
    
    global p1 p2 q1 q2 v
%     p1=1;p2=1;
%     q1=2;q2=3;
%     v=0.1;
    x1d=p1+v*(q1-p1)*t;
    x3d=atan2(q2-p2,q1-p1);
    x2d=p2+v*(q2-p2)*t;
    x1d_dot=v*(q1-p1);
    x3d_dot=0;
    x1d_ddot=0;
    x3d_ddot=0;
    eval1=5;eval2=6;
    eval3=5;eval4=6;
    v1=x1d_ddot-(eval1+eval2)*(y(2)-x1d_dot)-(eval1*eval2)*(y(1)-x1d);
    v2=x3d_ddot-(eval3+eval4)*(y(4)-x3d_dot)-(eval3*eval4)*(y(3)-x3d);
    
    % Control inputs 
    tau1=(a1*x(4)*cos(x(3))-a2*x(5)*sin(x(3))+(1-b2)*x(4)*x(6)*sin(x(3))...
        +(1-b1)*x(5)*x(6)*cos(x(3))+v1)/cos(x(3));
    tau2=a3*x(6)+b3*x(4)*x(5)+v2/c;
    
    % Dynamics
    dxdt=[x(4)*cos(x(3))-x(5)*sin(x(3));
        x(4)*sin(x(3))+x(5)*cos(x(3));
        c*x(6);
        -a1*x(4)+b1*x(5)*x(6)+tau1;
        -a2*x(5)-b2*x(4)*x(6);
        -a3*x(6)-b3*x(4)*x(5)+tau2];
    
end