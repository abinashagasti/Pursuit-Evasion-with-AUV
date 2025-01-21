classdef AUV_agent < handle

    properties
        
        % dynamics parameters based on system physical properties
        p11
        p12
        p21
        p22
        p31
        p32
        m1
        m2
        m3
        % coefficients of errors in sliding surfaces
        x_evals % [xe_dot xe integral xe]
        y_evals % [ye_dot ye integral ye]
        init_state % initial state
        % [x y psi u v r xd xd_dot xd_ddot xd_dddot yd yd_dot yd_ddot]
        T % total time
        pause_time % pause_time
        xd_ddddot % fourth derivative of desired x trajectory
        yd_dddot % third derivative of desired y trajectory
    end

    methods

        function agent = AUV_agent(model_parameters, initial_state, T, pause_time)
            agent.m1 = 1/(model_parameters.m-model_parameters.Xud);
            agent.m2 = 1/(model_parameters.m-model_parameters.Yvd);
            agent.m3 = 1/(model_parameters.Iz-model_parameters.Nrd);
            agent.p11 = agent.m1*model_parameters.Xu;
            agent.p12 = agent.m1*(model_parameters.m-model_parameters.Yvd);
            agent.p21 = agent.m2*model_parameters.Yv;
            agent.p22 = agent.m2*(model_parameters.Xud-model_parameters.m);
            agent.p31 = agent.m3*model_parameters.Nr;
            agent.p32 = agent.m3*(model_parameters.Yvd-model_parameters.Xud);
            agent.init_state = initial_state;
            agent.T = T;
            agent.pause_time = pause_time;
        end

        function set_trajectory(agent, xd_ddddot, yd_dddot)
            agent.xd_ddddot = xd_ddddot;
            agent.yd_dddot = yd_dddot;
        end

        function value = fdx(agent, state)
            value = -(state(10)+agent.x_evals(1)*state(9)+...
                agent.x_evals(2)*state(8)+agent.x_evals(3)*state(7));
        end

        function value = fdy(agent, t, state)
            value = -(agent.yd_dddot(t)+agent.y_evals(1)*state(13)+...
                agent.y_evals(2)*state(12)+agent.y_evals(3)*state(11));
        end

        function value = fdxd(agent, t, state)
            value = -(agent.xd_ddddot(t)+agent.x_evals(1)*state(10)+...
                agent.x_evals(2)*state(9)+agent.x_evals(3)*state(8));
        end

        function value = fs1(agent, state)
            value = agent.x_evals(1)*(agent.p11*state(4)+(agent.p12-1)*state(5)*state(6))*cos(state(3))+...
                agent.x_evals(2)*state(4)*cos(state(3))-agent.x_evals(1)*(agent.p21*state(5)+...
                (agent.p22+1)*state(4)*state(6))*sin(state(3))-agent.x_evals(2)*state(5)*sin(state(3))+...
                agent.x_evals(3)*state(1)+agent.fdx(state);
        end

        function value = fs2(agent, t, state)
            value = (agent.p11*(agent.p11*state(4)+agent.p12*state(5)*state(6))+...
                agent.p12*(agent.p21*state(5)*state(6)+agent.p22*state(4)*state(6)^2+...
                agent.p31*state(5)*state(6)+agent.p32*state(4)*state(5)^2)-...
                2*state(6)*(agent.p21*state(5)+agent.p22*state(4)*state(6))-...
                state(5)*(agent.p31*state(6)+agent.p32*state(4)*state(5))+...
                agent.y_evals(1)*(agent.p11*state(4)+(agent.p12-1)*state(5)*state(6))+...
                agent.y_evals(2)*state(4)-state(4)*state(6)^2)*sin(state(3))+...
                (agent.p21*(agent.p21*state(5)+agent.p22*state(4)*state(6))+...
                agent.p22*(agent.p11*state(4)*state(6)+agent.p12*state(5)*state(6)^2+...
                agent.p31*state(4)*state(6)+agent.p32*state(4)^2*state(5))+...
                2*state(6)*(agent.p11*state(4)+agent.p12*state(5)*state(6))+...
                state(4)*(agent.p31*state(6)+agent.p32*state(4)*state(5))+...
                agent.y_evals(1)*(agent.p21*state(5)+(agent.p22+1)*state(4)*state(6))+...
                agent.y_evals(2)*state(5)-state(5)*state(6)^2)*cos(state(3))+...
                agent.y_evals(3)*state(2)+agent.fdy(t,state);
        end

        function value = ftu(agent, state, tauu)
            value = agent.m1*tauu*((agent.p11+agent.y_evals(1))*sin(state(3))+ ...
                (agent.p22+2)*state(6)*cos(state(3)));
        end

        function value = ftu1d(agent, state, tauu)
            value = agent.m1*tauu*((agent.x_evals(1)*agent.p11+agent.x_evals(2))* ...
                cos(state(3))-agent.x_evals(1)*(agent.p22+1)*state(6)*sin(state(3)));
        end

        function value = fs11d(agent, t, state)
            value = ((agent.x_evals(1)*agent.p11+agent.x_evals(2))*(agent.p11*state(4) ...
                +agent.p12*state(5)*state(6))+agent.x_evals(1)*(agent.p12-1) ...
                *state(6)*(agent.p21*state(5)+agent.p22*state(4)*state(6))+ ...
                agent.x_evals(1)*(agent.p12-1)*state(5)*(agent.p31*state(6)+ ...
                agent.p32*state(4)*state(5)))*cos(state(3))- ...
                ((agent.x_evals(1)*agent.p11+agent.x_evals(2))*state(4)+ ...
                agent.x_evals(1)*(agent.p12-1)*state(5)*state(6))*state(6)*sin(state(3))- ...
                ((agent.x_evals(1)*agent.p21+agent.x_evals(2))*(agent.p21*state(5)+ ...
                agent.p22*state(4)*state(6))+agent.x_evals(1)*(agent.p22+1)*state(6)* ...
                (agent.p11*state(4)+agent.p12*state(5)*state(6))+agent.x_evals(1)*(agent.p22+1)* ...
                state(4)*(agent.p31*state(6)+agent.p32*state(4)*state(5)))*sin(state(3))- ...
                ((agent.x_evals(1)*agent.p21+agent.x_evals(2))*state(5)+ ...
                agent.x_evals(1)*(agent.p22+1)*state(4)*state(6))*state(6)*cos(state(3))+ ...
                agent.x_evals(3)*(state(4)*cos(state(3))-state(5)*sin(state(3)))+agent.fdxd(t,state);
            % value = -agent.fs1(state)*((1-agent.p22)*state(6)*tan(state(3))+ ...
            %     (agent.p11+agent.x_evals(2)/agent.x_evals(1)))+(agent.p21*state(5)+ ...
            %     agent.p22*state(4)*state(6))*((-agent.x_evals(1)*agent.p21- ...
            %     agent.x_evals(2))*sin(state(3))+(agent.x_evals(1)*agent.p11+ ...
            %     agent.x_evals(2))*cos(state(3)))+(agent.p11*state(4)+ ...
            %     agent.p12*state(5)*state(6))*((1-agent.p22)*agent.x_evals(1)* ...
            %     state(6)*sin(state(3))+(agent.x_evals(1)*agent.p11+agent.x_evals(2))* ...
            %     cos(state(3)))+(agent.p31*state(6)+agent.p32*state(4)*state(5))* ...
            %     ((1-agent.p22)*agent.x_evals(1)*state(4)*sin(state(3))+ ...
            %     (agent.p12-1)*agent.x_evals(1)*state(5)*cos(state(3)))- ...
            %     (agent.x_evals(1)*(agent.p11*state(4)*state(6)+(agent.p12-1)* ...
            %     state(5)*state(6)^2)+agent.x_evals(2)*state(4)*state(6)-agent.x_evals(3)*state(5)) ...
            %     *sin(state(3))-(agent.x_evals(1)*(agent.p21*state(5)*state(6)+(agent.p22-1)* ...
            %     state(4)*state(6)^2)+agent.x_evals(2)*state(5)*state(6)+agent.x_evals(3)*state(4)) ...
            %     *cos(state(3))+agent.fdxd(t,state);
            
        end

        function velocity = dynamics(agent, t, state)
            velocity = zeros(length(agent.init_state),1);
            tauu = -agent.fs1(state)/(agent.m1*agent.x_evals(1)*cos(state(3)));
            taur = ((-agent.fs11d(t,state)-agent.ftu1d(state,tauu)-agent.fs1(state)* ...
                state(6)*tan(state(3)))*tan(state(3)))/(agent.m3*agent.x_evals(1)* ...
                ((agent.p12-1)*state(5)-(agent.p22+1)*state(4))*cos(state(3))- ...
                ((agent.p22+1)*state(6)+(agent.p12-1)*state(5))*sin(state(3)));
            % taur = 0;
            velocity(1:6) = [state(4)*cos(state(3))-state(5)*sin(state(3));
                state(4)*sin(state(3))+state(5)*cos(state(3));
                state(6);
                agent.p11*state(4)+agent.p12*state(5)*state(6)+agent.m1*tauu;
                agent.p21*state(5)+agent.p22*state(4)*state(6)
                agent.p31*state(6)+agent.p32*state(4)*state(5)+agent.m3*taur;];
            velocity(7:10)=[state(8);state(9);state(10);agent.xd_ddddot(t)];
            velocity(11:13)=[state(12);state(13);agent.yd_dddot(t)];
        end

        function plot_trajectories(agent,t,state)
            figure
            plot(t,-state(:,1)+state(:,7))
            hold on
            plot(t,-state(:,2)+state(:,11))
            plot(t,-state(:,3)+atan2(state(:,2),state(:,1)))
            title('Error in States')
            xlabel('Time')
            ylabel('Error')
            % legend('x-error','y-error')
            legend('x-error','y-error','\psi-error')
            
            % figure
            % plot(t,state(:,7))
            % hold on
            % plot(t,tau2)
            % title('Control')
            % xlabel('Time')
            % ylabel('Value')
            % legend('tau1','tau2')
            
            figure;
            hold on;
            
            % Initialize the plot
            hTrajectory = plot(state(1,1), state(1,2), 'b'); % Trajectory line
            hPoint = plot(state(1,1), state(1,2), 'ro', 'MarkerFaceColor', 'r'); % Moving point
            
            % dTrajectory = plot(state(1,8), state(1,11), 'b'); % Trajectory line
            % dPoint = plot(state(1,1), state(1,2), 'ro', 'MarkerFaceColor', 'r'); % Moving point
            
            % Set axis limits
            xlim([min(state(:,1)), max(state(:,1))]);
            ylim([min(state(:,2)), max(state(:,2))]);
            
            % Animation loop
            for k = 2:length(t)
                % Update the trajectory
                set(hTrajectory, 'XData', state(1:k,1), 'YData', state(1:k,2));
                
                % Update the moving point
                set(hPoint, 'XData', state(k,1), 'YData', state(k,2));
                
                % Pause for a short duration to create the animation effect
                pause(agent.pause_time);
            end
            
            hold off;

        end

        function [t,state] = trajectory_computation(agent, xd_ddddot, yd_dddot, x_evals, y_evals, plot)
            agent.set_trajectory(xd_ddddot, yd_dddot);
            agent.x_evals = x_evals;
            agent.y_evals = y_evals;
            [t,state] = ode45(@agent.dynamics, [0,agent.T], agent.init_state);
            if plot
                agent.plot_trajectories(t,state);
            end
        end

    end


end