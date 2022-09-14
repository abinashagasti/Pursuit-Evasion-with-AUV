%% 
clear
clc
close all

%
p1=1;p2=1;
x0=[p1;p2;pi/4;0;0;0]; % Starting point
q1=20;q2=30; % Opponent location
v=1;ve=0.5;
param=[q1;q2;v;ve];
[t,x,xe,ye]=rk4(0.01,[0,6],x0,param);

%%
plot(x(1,:),x(2,:))
hold on
plot(xe,ye)
xlabel('x')
ylabel('y')
title('Actual system trajectory on xy plane')

% psid=atan2(q2-p2,q1-p1);
% figure
% plot(t,x(1,:)-(p1+v*cos(psid)*t))
% title('Error in x trajectory')
% % legend('x(t)','x_d(t)')
% figure
% plot(t,x(2,:)-(p2+v*sin(psid)*t))
% title('Error in y Trajectory')
% % legend('y(t)','y_d(t)')
% figure
% plot(t,x(3,:)-psid*ones(size(t)))
% title('Error in \psi Trajectory')
% % legend('\psi(t)','\psi_d(t)')