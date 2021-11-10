
classdef utilityFunctions
    methods (Static)
       
        %% 
        function A_Ni = change_system_representation(Ai,Fi,Ci,Agraph)
            M = size(Ai,2); % number of subsystems    
            G = digraph(Agraph);
            A_Ni = cell(1,M);
            for i=1:M
                out_neighbors = sort([i,successors(G, i)]);
                Acell = cell(1,length(out_neighbors));
                Acell{i} = Ai{i} - Fi{i}*sum(Agraph(i,:))*Ci{i};
                out_neighbors(i) = [];
                for j=1:length(out_neighbors)
                   Acell{out_neighbors(j)} = Fi{i}*Agraph(i, out_neighbors(j))*...
                                             Ci{out_neighbors(j)};
                end
                A_Ni{i} = cell2mat(Acell);    
            end
        end
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
        function plot_states_coupled_osci(X,U, config, control_type, param)  
            if config == "GENERAL"
                states = cell2mat(X); % extract position/velocity at each timesteps
                k = param.number_subsystem * param.size_subsystem;
                position{1} = states(1:k:end); position{2} = states(3:k:end);
                velocity{1} = states(2:k:end); velocity{2} = states(4:k:end); 
                controller = cell2mat(U');%first row u1, second row u2, column are timesteps
            elseif config == "DISTRIBUTED"
                k = param.size_subsystem;
                states = cell2mat(X);
                position{1} = states(1:k:end,1);
                position{2} = states(1:k:end,2);
                velocity{1} = states(2:k:end,1);
                velocity{2} = states(2:k:end,2);
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
        
        %%
        function plot_DGU_system(X,U, config, control_type, param)  
             M = param.number_subsystem;
             lgd = cell(1,M);
             voltage = cell(1,M); current= cell(1,M); integrator = cell(1,M);
             if config == "GENERAL"
                states = cell2mat(X); % extract position/velocity at each timesteps
                k = M * param.size_subsystem;
                j = 1;
                for i = 1:M
                    voltage{i} = states(j:k:end);
                    current{i} = states(j+1:k:end);
                    integrator{i} = states(j+2:k:end);
                    j = j+param.size_subsystem;
                    lgd{i} = sprintf("DGU %d", i);
                end
                controller = cell2mat(U');%first row u1, second row u2, column are timesteps
            elseif config == "DISTRIBUTED"
                k = param.size_subsystem;
                states = cell2mat(X);
                for i = 1:M
                    voltage{i} = states(1:k:end,i);
                    current{i} = states(2:k:end,i);
                    integrator{i} = states(3:k:end,i);
                    lgd{i} = sprintf("DGU %d", i);
                end
                controller = cell2mat(U)'; %first row u1, second row u2, column are timesteps
            else
                error("not implemented configuration in plot states");
            end

            figure()
            sgtitle(control_type);
            subplot(2,1,1)
            title('Voltages');
            
            hold on
            for i = 1:M
                plot(voltage{i});
            end
            legend(string(lgd));
            grid on
            hold off
            subplot(2,1,2)
            title('Converter Currents');
            hold on
            for i = 1:M
                plot(current{i});
            end
            legend(string(lgd));
            grid on
            hold off

            figure()
            title("Controller  " + control_type)
            hold on
            for i = 1:M
                plot(controller(i,:));
            end
            legend(string(lgd));
            grid on
            hold off
        end
    end % end methods
end % end class