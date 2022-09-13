function dxdt = demo_dynamics(t,x)
% This is a demo dynamics to test out the custom ODE solver

    dxdt=[0 1;-2 -3]*x;

end

