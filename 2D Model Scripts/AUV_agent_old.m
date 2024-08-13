classdef AUV_agent_old < handle
    properties
        model_parameters % parameters of given 2D model
        t % time of the simulation
        T % total simulation time
        del % del t interval
        % effective parameters of dynamic model
        a1
        a2
        a3
        b1
        b2
        b3
        % states of system
        x
        y
        psi
        u
        v
        r
        % controls of system
        tau1
        tau2
        tau1dot
        % stored variables
        variables
    end

    methods

        function agent = AUV_agent_old(model_parameters, initial_state, T, del)
            agent.model_parameters = model_parameters;
            agent.a1 = model_parameters.Xu/(model_parameters.m-model_parameters.Xud);
            agent.a2 = model_parameters.Yv/(model_parameters.m-model_parameters.Yvd);
            agent.a3 = model_parameters.Nr/(model_parameters.Iz-model_parameters.Nrd);
            agent.b1 = (model_parameters.m-model_parameters.Yvd)/(model_parameters.m-model_parameters.Xud);
            agent.b2 = (model_parameters.Xud-model_parameters.m)/(model_parameters.m-model_parameters.Yvd);
            agent.b3 = (model_parameters.Yvd-model_parameters.Xud)/(model_parameters.Iz-model_parameters.Nrd);
            agent.x = initial_state(1);
            agent.y = initial_state(2);
            agent.psi = initial_state(3);
            agent.u = initial_state(4);
            agent.v = initial_state(5);
            agent.r = initial_state(6);
            agent.tau1 = 0;
            agent.tau2 = 0;
            agent.tau1dot = 0;
            agent.t = 0;
            agent.T = T;
            agent.del = del;
            agent.variables = [agent.x;agent.y;agent.psi;agent.u;agent.v;agent.r;agent.tau1;agent.tau2];
        end


        function velocity = dynamics(agent, controller, xd, yd)
            % Define the dynamics of the system
            if strcmpi(controller,"feedback linearization")
                agent.feedback_linearization_controller(xd,yd);
            else
                error("The controller "+controller+" does not exist!")
            end
            velocity = [agent.u*cos(agent.psi)-agent.v*sin(agent.psi);
                        agent.u*sin(agent.psi)+agent.v*cos(agent.psi);
                        agent.r;
                        agent.a1*agent.u+agent.b1*agent.v*agent.r+agent.tau1;
                        agent.a2*agent.v+agent.b2*agent.u*agent.r;
                        agent.a3*agent.r+agent.b3*agent.u*agent.v+agent.tau2;
                        agent.tau1dot];
        end

        function [xddot,xdddot,xddddot,yddot,ydddot,yddddot] = get_desired_traj_stats(agent, xd, yd)
            tdur = linspace(0,agent.T,agent.T/agent.del);
            xddot = diff(xd)./diff(tdur);
            yddot = diff(yd)./diff(tdur);
            xdddot = [0,diff(xddot)]./diff(tdur);
            ydddot = [0,diff(yddot)]./diff(tdur);
            xddddot = [0,diff(xdddot)]./diff(tdur);
            yddddot = [0,diff(ydddot)]./diff(tdur);
            xddot = [0,xddot];
            yddot = [0,yddot];
            xdddot = [0,xdddot];
            ydddot = [0,ydddot];
            xddddot = [0,xddddot];
            yddddot = [0,yddddot];
        end

        function feedback_linearization_controller(agent, xd, yd)

            f1 = agent.a1^2*agent.u+((agent.a1+agent.a2+agent.a3)*agent.b1-2*agent.a2-agent.a3)*agent.v*agent.r+...
                (agent.b1*agent.b2-2*agent.b2-1)*agent.u*agent.r^2+(agent.b1-1)*agent.b3*agent.u*agent.v^2;
            f2 = agent.a2^2*agent.v+((agent.a1+agent.a2+agent.a3)*agent.b2+2*agent.a1+agent.a3)*agent.u*agent.r+...
                (agent.b1*agent.b2+2*agent.b1-1)*agent.v*agent.r^2+(agent.b2+1)*agent.b3*agent.u^2*agent.v;

            [xddot,xdddot,xddddot,yddot,ydddot,yddddot] = agent.get_desired_traj_stats(xd, yd);

            if agent.t/agent.del==0
                x_dot = 0;
                x_ddot = 0;
                y_dot = 0;
                y_ddot = 0;
            elseif agent.t/agent.del==1
                x_dot = (agent.variables(1,end)-agent.variables(1,end-1))/agent.del;
                x_ddot = 0;
                y_dot = (agent.variables(2,end)-agent.variables(2,end-1))/agent.del;
                y_ddot = 0;
            else
                x_dot = (agent.variables(1,end)-agent.variables(1,end-1))/agent.del;
                x_dot_1 = (agent.variables(1,end-1)-agent.variables(1,end-2))/agent.del;
                x_ddot = (x_dot_1-x_dot)/agent.del;
                y_dot = (agent.variables(2,end)-agent.variables(2,end-1))/agent.del;
                y_dot_1 = (agent.variables(2,end-1)-agent.variables(2,end-2))/agent.del;
                y_ddot = (y_dot_1-y_dot)/agent.del;
            end

            c1=3;c2=3;c3=1;c4=3;c5=3;c6=1;
            
            ux = xddddot(1+round(agent.t/agent.del))+c1*(xdddot(1+round(agent.t/agent.del))-x_ddot)+c2*(xddot(1+round(agent.t/agent.del))-x_dot)+c3*(xd(1+round(agent.t/agent.del))-agent.x);
            uy = yddddot(1+round(agent.t/agent.del))+c4*(ydddot(1+round(agent.t/agent.del))-y_ddot)+c5*(yddot(1+round(agent.t/agent.del))-y_dot)+c6*(yd(1+round(agent.t/agent.del))-agent.y);

            agent.tau2 = (1/((agent.b2+1)*agent.u))*(-(agent.b2+2)*agent.r*agent.tau1-f2-ux*sin(agent.psi)+uy*cos(agent.psi));
            agent.tau1dot = -(agent.a1*agent.tau1+(agent.b1-1)*agent.v*agent.tau2)-f1+ux*cos(agent.psi)+uy*sin(agent.psi);

        end

        function trajectory_computation(agent, xd, yd)
            for i=1:(agent.T/agent.del)
                velocity = agent.dynamics("feedback linearization", xd, yd);
                new_state = [agent.x;agent.y;agent.psi;agent.u;agent.v;agent.r;agent.tau1]+velocity*agent.del;
                agent.x = new_state(1);
                agent.y = new_state(2);
                agent.psi = new_state(3);
                agent.u = new_state(4);
                agent.v = new_state(5);
                agent.r = new_state(6);
                agent.tau1 = new_state(7);
                agent.t = agent.t+agent.del;
                agent.variables = [agent.variables,[agent.x;agent.y;agent.psi;agent.u;agent.v;agent.r;agent.tau1;agent.tau2]];
            end
        end

    end
end
