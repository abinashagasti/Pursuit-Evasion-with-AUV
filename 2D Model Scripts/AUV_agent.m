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
        % upper and lower bounds for controls
        tau1_M
        tau1_m
        tau2_M
        tau2_m

        T % total time duration
        pause_time % for plotting trajectory

        % initial states of system
        init_state 
        % [x y psi u v r tau1 xd xd_dot xd_ddot yd yd_dot yd_ddot]
        % 13 dimensions

        % desired trajectory
        xd_dddot
        yd_dddot

        % desired eigenvalues
        x_evals
        y_evals

        % polynomials
        x_coeffs
        y_coeffs

    end

    methods

        function agent = AUV_agent(model_parameters, initial_state, T, pause_time)
            agent.model_parameters = model_parameters;
            agent.a1 = model_parameters.Xu/(model_parameters.m-model_parameters.Xud);
            agent.a2 = model_parameters.Yv/(model_parameters.m-model_parameters.Yvd);
            agent.a3 = model_parameters.Nr/(model_parameters.Iz-model_parameters.Nrd);
            agent.b1 = (model_parameters.m-model_parameters.Yvd)/(model_parameters.m-model_parameters.Xud);
            agent.b2 = (model_parameters.Xud-model_parameters.m)/(model_parameters.m-model_parameters.Yvd);
            agent.b3 = (model_parameters.Yvd-model_parameters.Xud)/(model_parameters.Iz-model_parameters.Nrd);
            agent.init_state = initial_state;
            agent.T = T;
            agent.pause_time = pause_time;
            agent.tau1_M = Inf;
            agent.tau1_m = -Inf;
            agent.tau2_M = Inf;
            agent.tau2_m = -Inf;
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
            agent.x_coeffs = poly(agent.x_evals(t));
            c1=agent.x_coeffs(2);c2=agent.x_coeffs(3);c3=agent.x_coeffs(4);
            x_ddot = (agent.a1*state(4)+agent.b1*state(5)*state(6)+state(7))*cos(state(3))-state(4)*state(6)*sin(state(3))+...
                -(agent.a2*state(5)+agent.b2*state(4)*state(6))*sin(state(3))-state(5)*state(6)*cos(state(3));
            x_dot = state(4)*cos(state(3))-state(5)*sin(state(3));
            value = agent.xd_dddot(t)+c1*(state(10)-x_ddot)+c2*(state(9)-x_dot)+c3*(state(8)-state(1));
        end

        function value = uy(agent, t, state)
            agent.y_coeffs = poly(agent.y_evals(t));
            c4=agent.y_coeffs(2);c5=agent.y_coeffs(3);c6=agent.y_coeffs(4);
            y_ddot = (agent.a1*state(4)+agent.b1*state(5)*state(6)+state(7))*sin(state(3))+state(4)*state(6)*cos(state(3))+...
                +(agent.a2*state(5)+agent.b2*state(4)*state(6))*cos(state(3))-state(5)*state(6)*sin(state(3));
            y_dot = state(4)*sin(state(3))+state(5)*cos(state(3));
            value = agent.yd_dddot(t)+c4*(state(13)-y_ddot)+c5*(state(12)-y_dot)+c6*(state(11)-state(2));
        end

        function velocity = dynamics(agent, t, state)
            velocity = zeros(length(agent.init_state),1);
            tau2 = (1/((agent.b2+1)*state(4)))*(-(agent.b2+2)*state(6)*state(7)-agent.f2(state)-agent.ux(t,state)*sin(state(3))+agent.uy(t,state)*cos(state(3)));
            if tau2>agent.tau2_M
                tau2=agent.tau2_M;
            elseif tau2<agent.tau2_m
                tau2=agent.tau2_m;
            end
            velocity(1:6) = [state(4)*cos(state(3))-state(5)*sin(state(3));
                state(4)*sin(state(3))+state(5)*cos(state(3));
                state(6);
                agent.a1*state(4)+agent.b1*state(5)*state(6)+state(7);
                agent.a2*state(5)+agent.b2*state(4)*state(6)
                agent.a3*state(6)+agent.b3*state(4)*state(5)+tau2;];
                velocity(7)=-(agent.a1*state(7)+(agent.b1-1)*state(5)*tau2)-agent.f1(state)+agent.ux(t,state)*cos(state(3))+agent.uy(t,state)*sin(state(3));
                if state(7)>agent.tau1_M && velocity(7)>0
                    velocity(7)=0;
                elseif state(7)<agent.tau1_m && velocity(7)<0
                    velocity(7)=0;
                end
                velocity(8:10)=[state(9);state(10);agent.xd_dddot(t)];
                velocity(11:13)=[state(12);state(13);agent.yd_dddot(t)];
        end

        function tau2 = plot_trajectories(agent,t,state)
            figure
            plot(t,-state(:,1)+state(:,8))
            hold on
            plot(t,-state(:,2)+state(:,11))
            plot(t,-mod(state(:,3),2*pi)+mod(atan2(state(:,2),state(:,1)),2*pi))
            title('Error in States')
            xlabel('Time')
            ylabel('Error')
            % legend('x-error','y-error')
            legend('x-error','y-error','\psi-error')
            
            tau2 = zeros(length(t),1);
            for k=1:length(t)
                tau2(k) = (1/((agent.b2+1)*state(k,4)))*(-(agent.b2+2)*state(k,6)*state(k,7)-agent.f2(state(k,:))-agent.ux(t(k),state(k,:))*sin(state(k,3))+agent.uy(t(k),state(k,:))*cos(state(k,3)));
                if tau2(k)>agent.tau2_M
                    tau2(k)=agent.tau2_M;
                elseif tau2(k)<agent.tau2_m
                    tau2(k)=agent.tau2_m;
                end
            end
            figure
            plot(t,state(:,7))
            hold on
            plot(t,tau2)
            title('Control')
            xlabel('Time')
            ylabel('Value')
            legend('tau1','tau2')
            
            figure;
            % hold on;
            
            % Initialize the plot
            % hTrajectory = plot(state(1,1), state(1,2), 'b'); % Trajectory line
            % hPoint = plot(state(1,1), state(1,2), 'ro', 'MarkerFaceColor', 'r'); % Moving point
            % 
            % % dTrajectory = plot(state(1,8), state(1,11), 'b'); % Trajectory line
            % % dPoint = plot(state(1,1), state(1,2), 'ro', 'MarkerFaceColor', 'r'); % Moving point
            % 
            % % Set axis limits
            % xlim([min(state(:,1)), max(state(:,1))]);
            % ylim([min(state(:,2)), max(state(:,2))]);
            % 
            % % Animation loop
            % for k = 2:length(t)
            %     % Update the trajectory
            %     set(hTrajectory, 'XData', state(1:k,1), 'YData', state(1:k,2));
            % 
            %     % Update the moving point
            %     set(hPoint, 'XData', state(k,1), 'YData', state(k,2));
            % 
            %     % Pause for a short duration to create the animation effect
            %     pause(agent.pause_time);
            % end

            plot(state(:,1),state(:,2))
            hold on
            plot(state(:,8),state(:,11))
            title('AUV Trajectory')
            xlabel('x-axis')
            ylabel('y-axis')
            legend('State Trajectory','Desired Trajectory')
            
            hold off;

        end

        function [t,state,tau2] = trajectory_computation(agent, xd_dddot, yd_dddot, x_evals, y_evals, plot)
            agent.set_trajectory(xd_dddot, yd_dddot);
            agent.x_evals = x_evals;
            agent.y_evals = y_evals;
            [t,state] = ode45(@agent.dynamics, [0,agent.T], agent.init_state);
            if plot
                tau2 = agent.plot_trajectories(t,state);
            end
        end

    end


end