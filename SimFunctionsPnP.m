classdef SimFunctionsPnP
    methods (Static)
        %% ------------------ COMPUTE PASSIVE GAINS ---------------------------%
        function obj = setPassiveControllers(obj)
            for i= obj.activeDGU
                [obj.Ki{i}, Di, obj.Pi{i}, Gamma_i] = controller_passivity(...
                                             obj.Ai{i}, obj.Bi{i},...
                                             obj.Ci{i}, obj.Fi{i}, ...
                                             obj.L_tilde, obj.global_sysd.C,i);
                sprintf("K%d", i)
                disp(obj.Ki{i});
                sprintf("P%d", i)
                disp(obj.Pi{i});
                fprintf("minimum eigenvalue of dissipation rate %.2d \n", ...
                        min(eig(Gamma_i)));
               
            end
            for i= obj.activeDGU
                neighbors_i = sort([i; neighbors(obj.NetGraph, i)]); 
                K_block = blkdiag(obj.Ki{neighbors_i}); %block matrix of all neighbors K
                obj.K_Ni{i} = K_block(neighbors_i==i,:); % extract only row corresponding to subsystem i
                %[Q_Ni, Ri] = computeQi_Ri(obj, i);
            end
            
        end
        
        
        %% -------------SIMPLE SIMULATIONS IN DELTA FORMULATION----------------%
        % -------- SIMULATION OF DGU NETWORK FOR LOCAL CONTROLLERS Ki ---------%
        function [X,U] = sim_DGU_distributed(x0, length_sim, param, K)
            M = param.nb_subsystems; % number of subsystems
            X = cell(1,length_sim+1,1); % state at each timestep
            U = cell(1, length_sim); % control input at each timestep
            X{1} = horzcat(x0{:}); % columns are subsystem i

            if ~param.delta_config
                error("Switch to delta configuration to use Passivity to converge to ref");
            end

            dX = cell(length_sim+1,1); % state at each timestep
            dU = cell(length_sim,1); % control input at each timestep
            dX{1} = X{1} - horzcat(param.Xref{:});

            for n=1:length_sim
                 for i=1:M
                    dU{n}(:,i) = K{i}*dX{n}(:,i);
                    U{n} = dU{n} + param.Uref{i};
                    out_neighbors = sort([i; neighbors(param.NetGraph, i)]); % get neighbors
                    % create neighbor state vector comprising the state of subsystem i
                    % and the state of it's neighbors (concatenated)
                    x_Ni = reshape(dX{n}(:,out_neighbors),[],1); 
                    % Apply first input to the system and get next state
                    dX{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*dU{n}(:,i);
                    X{n+1}(:,i) = dX{n+1}(:,i) + param.Xref{i};
                 end
            end  
        end
        
        % -------- SIMULATION OF DGU NETWORK FOR GLOBAL CONTROLLER K ---------%
        function [X,U] = sim_global_DGU(x0, length_sim, param, K)
            Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
            X = cell(length_sim+1,1); % state at each timestep
            U = cell(length_sim,1); % control input at each timestep
            X{1} = vertcat(x0{:}); % global initial state, vector Mx1
            if ~param.delta_config
                error("Switch to delta configuration to use Passivity to converge to ref");
            end
            dX = cell(length_sim+1,1); % state at each timestep
            dU = cell(length_sim,1); % control input at each timestep
            dX{1} = X{1} - vertcat(param.Xref{:});
            for n=1:length_sim
                dU{n} = K*dX{n};
                U{n} = dU{n} + vertcat(param.Uref{:});
                dX{n+1}= Ad*dX{n}+ Bd*dU{n};
                X{n+1} = dX{n+1} + vertcat(param.Xref{:});
            end
        end
        %% ------------------------ P&P OPERATIONS------------------------------
        %-------------------------- REDESIGN PHASE-----------------------------%       
        function obj = redesignPhase(obj, oldGraph, idxDGU, procedure)
            if procedure == "add"
                out_neighbors = sort([idxDGU; neighbors(oldGraph, idxDGU)]);  
            elseif procedure == "delete"
                out_neighbors = sort(neighbors(oldGraph, idxDGU));
            else
                error("procedure not well defined: add or delete");
            end
             for k = 1:length(out_neighbors)
                 i = out_neighbors(k);
                 [obj.Ki{i}, Di, obj.Pi{i}, Gamma_i] = controller_passivity(...
                                                 obj.Ai{i}, obj.Bi{i},...
                                                 obj.Ci{i}, obj.Fi{i}, ...
                                                 obj.L_tilde, obj.global_sysd.C, i);
                sprintf("New passive controller gain of system %d", i)
                disp(obj.Ki{i});
                sprintf("P%d", i)
                disp(obj.Pi{i});
                fprintf("minimum eigenvalue of dissipation rate %.2d \n", ...
                        min(eig(Gamma_i)));
             end 
             for i= obj.activeDGU % recompute neighbors passive gains
                out_neighbors = sort([i; neighbors(obj.NetGraph, i)]); 
                K_block = blkdiag(obj.Ki{out_neighbors}); %block matrix of all neighbors K
                obj.K_Ni{i} = K_block(out_neighbors==i,:); % extract only row corresponding to subsystem i
             end
        end
        
    
        %-----------------------TRANSITION PHASE-------------------------------%
        function [X, U, lenSim, xs, us, alpha] = transitionPhase(x0, paramBefore,...
                                                         paramAfter, Qi, Ri, target, ADMM)
            N = 5; %Horizon 
            % compute steady state
            if ADMM
                disp("Using ADMM !");
                [xs, us, alpha] = transition_compute_ss_admm(horzcat(x0{:}), N, paramBefore,...
                                    paramAfter, target);
            elseif ~ADMM
                disp("Not using ADMM !");
                [xs, us, alpha] = transition_compute_ss(horzcat(x0{:}), N, paramBefore,...
                                    paramAfter, target);
            else
                disp("Error: Choose ADMM to be true or false");
                
            end

            % MPC REGULATION to steady state
            X{1} = horzcat(x0{:});
            clear regulation2ss
            n = 1;
            while any(abs(X{n}(1,paramBefore.activeDGU) - xs(1,paramBefore.activeDGU)) > 1e-1) || ...
                  any(abs(X{n}(2,paramBefore.activeDGU) - xs(2,paramBefore.activeDGU)) > 5e-2)    
                    % control input is of size nu x M 
                    if ADMM
                        U{n} = regulation2ss_admm(X{n}, N, paramBefore, xs, us, Qi, Ri); % get first control input
                    else
                        U{n} = regulation2ss(X{n}, N, paramBefore, xs, us, Qi, Ri); % get first control input
                    end
                    if isnan(U{n})
                        error("Input to apply to controller is Nan at iteration %d",n);
                    end
                    for i=paramBefore.activeDGU % Only modify the activated DGU
                        neighbors_i = [i; neighbors(paramBefore.NetGraph, i)]; % get neighbors
                        neighbors_i = sort(neighbors_i); % sorted neighbor list
                        % create neighbor state vector comprising the state of subsystem i
                        % and the state of it's neighbors (concatenated)
                        x_Ni = reshape(X{n}(:,neighbors_i),[],1); 
                        X{n+1}(:, i) = paramBefore.A_Ni{i}*x_Ni + paramBefore.Bi{i}*U{n}(:,i);
                    end
                    fprintf("Regulation to steady-state: iteration %d \n", n);
                    n = n+1;
            end  
            lenSim = n;
            if lenSim == 1 % if already at steady state
                U{n} = us;
            end
            for i=1:paramBefore.nb_subsystems
                if ~any(i == paramBefore.activeDGU(:))
                    X{1}(:,i) = [NaN;NaN];
                    for n = 1:lenSim
                        X{n+1}(:,i) = [NaN;NaN];
                        U{n}(:,i) = NaN;
                    end
                end
            end                                   
        end
        
        %% ------------------------- MPC CONTROLLERS--------------------------%%
        %----- TRACKING MPC WITH RECONFIGURABLE TERMINAL INGREDIENTS ----------%
        function [X,U] = mpc_DGU_tracking(controller, x0, length_sim,param, Q_Ni,...
                                          Ri)
                                      
            X = cell(1,length_sim+1,1); % state at each timestep
            U = cell(1, length_sim); % control input at each timestep
            %X{1} = horzcat(x0{:}); % initial state columns are subsystem i
            % initialize all steps and states, otherwise yalmip returns NaN in
            % non-activated DGU, resulting in NaN in the states which are given
            % again to yalmip which does not support NaN as x0
            X(:) = {horzcat(x0{:})}; 
            controllerType = functions(controller);
            fprintf("--INFO: Using controller %s -- \n ", controllerType.function);
            clear (controllerType.file)
             
            N = 5; % Horizon
       
            for n = 1:length_sim % loop over all subsystem
                % control input is of size nu x M 
                U{n} = controller(X{n}, Q_Ni, Ri, N, param); % get first control input
                if isnan(U{n})
                    error("Input to apply to controller is Nan at iteration %d",n);
                end
                for i=param.activeDGU % Only modify the activated DGU
                    neighbors_i = [i; neighbors(param.NetGraph, i)]; % get neighbors
                    neighbors_i = sort(neighbors_i); % sorted neighbor list
                    % create neighbor state vector comprising the state of subsystem i
                    % and the state of it's neighbors (concatenated)
                    x_Ni = reshape(X{n}(:,neighbors_i),[],1); 
                    X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
                end
            end
            for i=1:param.nb_subsystems
                % replace the non-active DGU states by NaN (for practical reason
                % when states are plotted afterwards).
                if ~any(i == param.activeDGU(:)) 
                    X{1}(:,i) = [NaN;NaN];
                    for n = 1:length_sim 
                        X{n+1}(:,i) = [NaN;NaN];
                        U{n}(:,i) = NaN;
                    end
                end
            end
        end
        
        %----- MPC DELTA FORMULATION WITH OFFLINE TERMINAL INGREDIENTS --------%
        function [X,U] = mpc_sim_DGU_delta(controller, x0, length_sim, param,...
                                          alpha, Q_Ni, Ri, Pi, Gamma_Ni)
                                      
            M = param.nb_subsystems; % number of subsystems
            dX = cell(1,length_sim+1); % state at each timestep
            dU = cell(1,length_sim); % control input at each timestep
            X = cell(1,length_sim+1); % state at each timestep
            U = cell(1,length_sim); % control input at each timestep
            X{1} = horzcat(x0{:}); % initial state columns are subsystem i    
            dX{1} = X{1} - horzcat(param.Xref{:});
            N = 10; % Horizon
            clear mpc_online
            for n = 1:length_sim % loop over all subsystem
                % control input is of size nu x M
                dU{n} = controller(dX{n}, alpha, Q_Ni, Ri, Pi, N, param); % get first control input
                U{n} = dU{n} + horzcat(param.Uref{:});
                if isnan(dU{n})
                    error("Input to apply to controller is Nan at iteration %d",n);
                end
                for i=1:M
                    out_neighbors = [i; neighbors(param.NetGraph, i)]; % get neighbors
                    out_neighbors = sort(out_neighbors); % sorted neighbor list
                    % create neighbor state vector comprising the state of subsystem i
                    % and the state of it's neighbors (concatenated)
                    x_Ni = reshape(dX{n}(:,out_neighbors),[],1); 
                    % Apply first input to the system and get next state
                    dX{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*dU{n}(:,i);
                    X{n+1}(:,i) = dX{n+1}(:,i) + param.Xref{i};
                    % update variable constraining terminal set of each subsystem
                    alpha(i) = alpha(i) + x_Ni'*Gamma_Ni{i}*x_Ni;
                end
            end                                                                             
        end
        
    end % end of static methods
end % class end