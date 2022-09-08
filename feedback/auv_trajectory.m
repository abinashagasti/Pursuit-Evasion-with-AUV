clear
clc
close all

%%
[t,x] = ode45(@auv_dynamics,[0,15],[1,1,pi/4,0,0,0]);
% plot(t,x(:,1),'-o',t,x(:,2),'-o')
% plot(t,x(:,1:2))
% hold on
plot(x(:,1),x(:,2),'b')
hold on
plot(3,3,'ro')
title('AUV Tracking Trajectories');
xlabel('x axis');
ylabel('y axis');

%%
figure
plot(t,x(:,1:6))
legend('x','y','\psi','u','v','r')