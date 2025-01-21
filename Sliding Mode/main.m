clear
clc
close all

%%
set(0,'DefaultFigureWindowStyle','docked')
% set(0,'DefaultFigureWindowStyle','normal')

model_parameters = struct("m",30.48,"Iz",3.45,"Xu",-8.8065,"Yv",-65.5457,"Nr",-6.7352,"Xud",-0.93,"Yvd",-35.5,"Nrd",-35.5);
%model_parameters = struct("p11",-0.2804,"p12",2.1006,"p21",-0.9934,"p22",-0.4761,"p31",-0.1729,"p32",-0.8875,"m1",0.0318,"m2",0.0257);

T = 15; % Total time
pause_time = 0.01;

% Define desired trajectories
xd_ddddot = @(t) 0; 
yd_dddot = @(t) 0;
% xd(t) = t; yd(t) = t;

% xd_dddot = @(t) sin(t); 
% yd_dddot = @(t) -cos(t);
% xd(t) = cos t; yd(t) = sin t;

initial_state = [0;0;pi/4;1;0;0; % states
    0;1;0;0; % xd(0); xd_dot(0); xd_ddot(0); xd_dddot(0)
    0;1;0]; % yd(0); yd_dot(0); yd_ddot(0)

agent = AUV_agent(model_parameters,initial_state,T,pause_time);


%%
x_evals = -1*ones(1,3);
y_evals = -1*ones(1,3);

[t,state] = agent.trajectory_computation(xd_ddddot,yd_dddot,x_evals,y_evals,true);
