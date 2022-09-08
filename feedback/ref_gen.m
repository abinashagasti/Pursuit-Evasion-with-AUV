% function [x,y,psi]=ref_gen(p1,p2,q1,q2,v)
% Generates a reference straight line trajectory 
    t=linspace(0,1);
    x=p1+v*(q1-p1)*t;
    y=p2+v*(q2-p2)*t;
    psi=atan2(y,x);

    plot(x,y)
    xlabel('x')
    ylabel('y')
    title('Desired Trajectory')
    figure
    plot(t,x)
    hold on
    plot(t,y)
    plot(t,psi)
    title('Desired Output States')
    legend('x','y','\psi')
% end