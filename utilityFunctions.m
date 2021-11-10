
classdef utilityFunctions
    methods (Static)
       
        %% Compute Finite Cost given state, input, Q and R
        function cost = compute_QR_cost(X,U,Q,R, config)
            cost = 0;
            if config == "GENERAL"
                for i=1:length(X)-1
                    cost = cost + X{i}'*Q*X{i} + U{i}'*R*U{i};
                end
            elseif config == "DISTRIBUTED"
                for i=1:length(X)-1
                    x = reshape(X{i},[],1); u = U{i}';
                    cost = cost + x'*Q*x + u'*R*u;
                end
            else
                error("wrong configuration, choose between general or distributed");
            end   
        end
        %% Plot states and control input for the coupled oscillator model
        function plot_states_coupled_osci(X,U, config, control_type)  
            if config == "GENERAL"
                states = cell2mat(X);
                position{1} = states(1:4:end); position{2} = states(3:4:end);
                velocity{1} = states(2:4:end); velocity{2} = states(4:4:end); 
                controller = cell2mat(U');%first row u1, second row u2, column are timesteps
            elseif config == "DISTRIBUTED"
                states = cell2mat(X);
                position{1} = states(1:2:end,1);
                position{2} = states(1:2:end,2);
                velocity{1} = states(2:2:end,1);
                velocity{2} = states(2:2:end,2);
                controller = cell2mat(U)'; %first row u1, second row u2, column are timesteps
            else
                error("not implemented configuration in plot states");
            end

            figure()
            sgtitle(control_type);
            subplot(2,1,1)
            title('Positions');
            hold on
            plot(position{1}, 'r-');
            plot(position{2}, 'b-');
            legend("mass 1", "mass 2");
            grid on
            hold off
            subplot(2,1,2)
            title('velocities');
            hold on
            plot(velocity{1}, 'r-');
            plot(velocity{2}, 'b-');
            legend("mass 1", "mass 2");
            grid on
            hold off

            figure()
            title("Controller  " + control_type)
            hold on
            plot(controller(1,:), 'r-');
            plot(controller(2,:), 'b-');
            legend("U1", "U2");
            grid on
            hold off
        end
        
    end % end methods
end % end class