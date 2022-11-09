clear
clc
close all

%%

global p1 p2 q1 q2 v

p1=0;p2=0;
q1=2;q2=3;
v=0.1;

[t,x] = ode45(@auv_dynamics,[0,20],[p1,p2,pi/4,0,0,0]);
% plot(t,x(:,1),'-o',t,x(:,2),'-o')
% plot(t,x(:,1:2))
% hold on
plot(x(:,1),x(:,2),'b')
hold on
plot(q1,q2,'ro')
title('AUV Tracking Trajectories');
xlabel('x axis');
ylabel('y axis');

%%
figure
plot(t,x(:,1:6))
legend('x','y','\psi','u','v','r')

%%
figure
plot(t,x(:,1)-(p1+v*(q1-p1)*t))
title('Error in x trajectory')
% legend('x(t)','x_d(t)')
figure
plot(t,x(:,2)-(p2+v*(q2-p2)*t))
title('Error in y Trajectory')
% legend('y(t)','y_d(t)')
figure
plot(t,x(:,3)-atan2(q2-p2,q1-p1)*ones(size(t)))
title('Error in \psi Trajectory')
% legend('\psi(t)','\psi_d(t)')

%%
figure
plot(x(:,1),x(:,2))
hold on
plot(q1+(v/2)*(q1-p1)*t,q2+(v/2)*(q2-p2)*t)



