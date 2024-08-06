clear
clc
close all

[t1,y1] = ode45(@(t,y1) 2*t+3+2*(t^2+3*t+5-y1),[0,10],0);

[t2,y2] = ode45(@(t,y2) 2*t+3+2*(t^2+3*t+5-y2)+10*rand,[0,10],0);

plot(t1,y1)
hold on
plot(t2,y2)
plot(t1,t1.^2+3*t1+5)
legend('Without noise','With noise','Desired trajectory')