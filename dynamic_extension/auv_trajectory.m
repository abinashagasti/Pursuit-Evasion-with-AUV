%% Computing the AUV dynamics 
% Tracking a straight line to pursue a holonomic evader

clear
clc
close all

global p1 p2 q1 q2 v

p1=0;p2=0;
q1=20;q2=30;
v=0.01;

% Run for pursuit evasion
[t,x] = ode45(@auv_dynamics,[0,200],[p1;p2;0;1;0;0;0]);
% plot(t,x(:,1),'-o',t,x(:,2),'-o')
% plot(t,x(:,1:2))
% hold on
plot(x(:,1),x(:,2),'b')
hold on
plot(q1,q2,'ro')
title('AUV Tracking Trajectories');
xlabel('x axis');
ylabel('y axis');

% Run for trajectory tracking
% [t,x] = ode45(@auv_dynamics,[0,200],[p1;p2;pi/4;1;0;0;0]);
% % plot(t,x(:,1),'-o',t,x(:,2),'-o')
% % plot(t,x(:,1:2))
% % hold on 
% plot(x(:,1),x(:,2),'b')
% hold on
% %plot(q1,q2,'ro')
% plot(t/10,sin(t/10),'r')
% title('AUV Tracking Trajectories');
% xlabel('x axis');
% ylabel('y axis');


%% Error in state trajectories

figure
plot(t,x(:,1)-(p1+v*(q1-p1)*t))
hold on
legend('Error in x')
% figure
plot(t,x(:,2)-(p2+v*(q2-p2)*t))
legend('Error in y')
% figure
plot(t,x(:,3)-atan2(q2-p2,q1-p1)*ones(size(t)))
legend('Error in x','Error in y','Error in \psi')
title('Error in state trajectories')
% legend('\psi(t)','\psi_d(t)')

%% Trajectory in 2D plane

figure
plot(x(:,1),x(:,2))
hold on
plot(q1+(v/2)*(q1-p1)*t,q2+(v/2)*(q2-p2)*t,'r')
plot(p1,p2,'bo')
plot(q1,q2,'ro')
title('Trajectories in 2D Plane')
legend('Pursuer','Evader')

%% Storing plot data

% writematrix(x(:,1),"/Users/abinashagasti/Work/IITM/Differential Games/Reports/AUV/xdata.dat")
% writematrix(x(:,2),"/Users/abinashagasti/Work/IITM/Differential Games/Reports/AUV/ydata.dat")
% writematrix(x(:,3),"/Users/abinashagasti/Work/IITM/Differential Games/Reports/AUV/psidata.dat")
% writematrix(x(:,1)-(p1+v*(q1-p1)*t),"/Users/abinashagasti/Work/IITM/Differential Games/Reports/AUV/xerror.dat")
% writematrix(x(:,2)-(p2+v*(q2-p2)*t),"/Users/abinashagasti/Work/IITM/Differential Games/Reports/AUV/yerror.dat")
% writematrix(x(:,3)-atan2(q2-p2,q1-p1)*ones(size(t)),"/Users/abinashagasti/Work/IITM/Differential Games/Reports/AUV/psierror.dat");
