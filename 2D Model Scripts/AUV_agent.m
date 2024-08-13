classdef AUV_agent < handle

    properties

        model_parameters % parameters of given 2D model
        % effective parameters of dynamic model
        a1
        a2
        a3
        b1
        b2
        b3

        T % total time duration

        % states of system
        state 
        % [x y psi u v r tau1 xd xd_dot xd_ddot yd yd_dot yd_ddot]
        % 13 dimensions

        % desired trajectory
        xd_dddot
        yd_dddot

    end

    methods

        function agent = AUV_agent(model_parameters, initial_state, T)
            agent.model_parameters = model_parameters;
            agent.a1 = model_parameters.Xu/(model_parameters.m-model_parameters.Xud);
            agent.a2 = model_parameters.Yv/(model_parameters.m-model_parameters.Yvd);
            agent.a3 = model_parameters.Nr/(model_parameters.Iz-model_parameters.Nrd);
            agent.b1 = (model_parameters.m-model_parameters.Yvd)/(model_parameters.m-model_parameters.Xud);
            agent.b2 = (model_parameters.Xud-model_parameters.m)/(model_parameters.m-model_parameters.Yvd);
            agent.b3 = (model_parameters.Yvd-model_parameters.Xud)/(model_parameters.Iz-model_parameters.Nrd);
            agent.state = initial_state;
            agent.T = T;
        end

        function set_trajectory(agent, xd_dddot, yd_dddot)
            agent.xd_dddot = xd_dddot;
            agent.yd_dddot = yd_dddot;
        end

        function value = f1(agent, state)
            value = agent.a1^2*state(4)+((agent.a1+agent.a2+agent.a3)*agent.b1-2*agent.a2-agent.a3)*state(5)*state(6)+...
                (agent.b1*agent.b2-2*agent.b2-1)*state(4)*state(6)^2+(agent.b1-1)*agent.b3*state(4)*state(5)^2;
        end

        function value = f2(agent,state)
            value = agent.a2^2*state(5)+((agent.a1+agent.a2+agent.a3)*agent.b2+2*agent.a1+agent.a3)*state(4)*state(6)+...
                (agent.b1*agent.b2+2*agent.b1-1)*state(5)*state(6)^2+(agent.b2+1)*agent.b3*state(4)^2*state(6);
        end

        function value = ux(agent, t, state)
            c1=3;c2=3;c3=1;
            x_ddot = (agent.a1*state(4)+agent.b1*state(5)*state(6)+state(7))*cos(state(3))-state(4)*state(6)*sin(state(3))+...
                -(agent.a2*state(5)+agent.b2*state(4)*state(6))*sin(state(3))-state(5)*state(6)*cos(state(3));
            x_dot = state(4)*cos(state(3))-state(5)*sin(state(3));
            value = agent.xd_dddot(t)+c1*(state(10)-x_ddot)+c2*(state(9)-x_dot)+c3*(state(8)-state(1));
        end

        function value = uy(agent, t, state)
            c4=3;c5=3;c6=1;
            y_ddot = (agent.a1*state(4)+agent.b1*state(5)*state(6)+state(7))*sin(state(3))+state(4)*state(6)*cos(state(3))+...
                +(agent.a2*state(5)+agent.b2*state(4)*state(6))*cos(state(3))-state(5)*state(6)*sin(state(3));
            y_dot = state(4)*sin(state(3))+state(5)*cos(state(3));
            value = agent.yd_dddot(t)+c4*(state(13)-y_ddot)+c5*(state(12)-y_dot)+c6*(state(11)-state(2));
        end

        function velocity = dynamics(agent, t, state)
            velocity = zeros(length(agent.state),1);
            tau2 = (1/((agent.b2+1)*state(4)))*(-(agent.b2+2)*state(6)*state(7)-agent.f2(state)-agent.ux(t,state)*sin(state(3))+agent.uy(t,state)*cos(state(3)));
            velocity(1:6) = [state(4)*cos(state(3))-state(5)*sin(state(3));
                state(4)*sin(state(3))+state(5)*cos(state(3));
                state(6);
                agent.a1*state(4)+agent.b1*state(5)*state(6)+state(7);
                agent.a2*state(5)+agent.b2*state(4)*state(6)
                agent.a3*state(6)+agent.b3*state(4)*state(5)+tau2;];
                velocity(7)=-(agent.a1*state(7)+(agent.b1-1)*state(5)*tau2)-agent.f1(state)+agent.ux(t,state)*cos(state(3))+agent.uy(t,state)*sin(state(3));
                velocity(8:10)=[state(9);state(10);agent.xd_dddot(t)];
                velocity(11:13)=[state(12);state(13);agent.yd_dddot(t)];
        end

        function [t,state] = trajectory_computation(agent, xd_dddot, yd_dddot)
            agent.set_trajectory(xd_dddot, yd_dddot);
            [t,state] = ode45(@agent.dynamics, [0,agent.T], agent.state);
        end

    end


end