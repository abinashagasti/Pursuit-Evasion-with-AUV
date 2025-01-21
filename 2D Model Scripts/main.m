clear
clc
close all

% Initialize model and desired trajectories

set(0,'DefaultFigureWindowStyle','docked')
% set(0,'DefaultFigureWindowStyle','normal')

model_parameters = struct("m",30.48,"Iz",3.45,"Xu",-8.8065,"Yv",-65.5457,"Nr",-6.7352,"Xud",-0.93,"Yvd",-35.5,"Nrd",-35.5);

T = 10*pi; % Total time
pause_time = 0.01;

% Define desired trajectories
xd_dddot = @(t) 0; 
yd_dddot = @(t) 0;
% xd(t) = t; yd(t) = t;

initial_state = [0;0;0;0.3;0;0; % states
    0; % tau1
    0;1;0; % xd(0); xd_dot(0); xd_ddot(0)
    0;1;0]; % yd(0); yd_dot(0); yd_ddot(0)

% xd_dddot = @(t) sin(t); 
% yd_dddot = @(t) -cos(t);
% % xd(t) = cos t; yd(t) = sin t;
% 
% initial_state = [1;0;pi/2;0.3;0;0; % states
%     0; % tau1
%     1;0;-1; % xd(0); xd_dot(0); xd_ddot(0)
%     0;1;0]; % yd(0); yd_dot(0); yd_ddot(0)

% Initialize AUV agent
% tau1_bound = 10;
% tau2_bound = -10;

agent = AUV_agent(model_parameters,initial_state,T,pause_time);

eval_bound = 10;
x_evals = @(t) -((t<=eval_bound).*t+(t>eval_bound).*eval_bound)*ones(1,3); 
y_evals = @(t) -((t<=eval_bound).*t+(t>eval_bound).*eval_bound)*ones(1,3);
% x_evals = @(t) -10*ones(1,3);
% y_evals = @(t) -10*ones(1,3);

[t,state,tau2] = agent.trajectory_computation(xd_dddot,yd_dddot,x_evals,y_evals,true);

 