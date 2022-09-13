function [t,x] = rk4(h,T,x0,param)
% This is a custom RK order 4 ODE solver that allows to tinker with the
% code better.
%     h=0.1;
    t=T(1):h:T(2);
    f=@auv_dynamics_adaptive;
    x=zeros(length(x0),length(t));
    x(:,1)=x0;
    param=[x0(1);x0(2);param(1:3)]; %Switch for point tracking
%     q1=param(1);q2=param(2);v=param(3);ve=param(4);
%     psid=atan2(q2-x0(2),q1-x0(1));
%     param=[x(1,1);x(2,1);q1;q2;v];
    for i=1:length(t)-1
%         xe=q1+ve*cos(psid)*t(i);
%         ye=q2+ve*sin(psid)*t(i);
%         if mod(i,20)==0
%            param=[x(1,i);x(2,i);xe;ye;v];
%         end
        k1=h*f(t(i),x(:,i),param);
        k2=h*f(t(i)+h/2,x(:,i)+k1/2,param);
        k3=h*f(t(i)+h/2,x(:,i)+k2/2,param);
        k4=h*f(t(i)+h,x(:,i)+k3,param);
        k=(k1+2*k2+2*k3+k4)/6;
        x(:,i+1)=x(:,i)+k;
    end

end

