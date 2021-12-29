classdef SimFunctionsPnP
    methods (Static)
        %% ------------------ COMPUTE PASSIVE GAINS ---------------------------%
        function [obj, lambda] = setPassiveControllers(obj)
            for i= obj.activeDGU
                [obj.Ki{i}, Di, obj.Pi{i}, Gamma_i] = controller_passivity(...
                                             obj.Ai{i}, obj.Bi{i},...
                                             obj.Ci{i}, obj.Fi{i}, ...
                                             obj.L_tilde, obj.global_sysd.C,i);
                sprintf("K%d", i)
                disp(obj.Ki{i});
                sprintf("P%d", i)
                disp(obj.Pi{i});
                lambda = min(eig(Gamma_i));
                fprintf("minimum eigenvalue of dissipation rate %.2d \n", ...
                        lambda);
               
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
        function [obj, Qi, Ri, Q_Ni] = redesignPhase(obj, oldGraph, idxDGU, procedure, QiOld, RiOld)
            if procedure == "add"
                out_neighbors = sort([idxDGU; neighbors(oldGraph, idxDGU)]);  
            elseif procedure == "delete"
                out_neighbors = sort(neighbors(oldGraph, idxDGU));
            else
                error("procedure not well defined: add or delete");
            end
            Qi = QiOld; Ri = RiOld;
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
             
             for k = 1:length(out_neighbors)
                % Compute associated Qi and RI matrices
                [Qi{i}, Ri{i}, decVariables{i}] = computeQi_Ri(obj, i);
                decVarSum = sum(cell2mat(decVariables),2);
                if all(decVarSum <= 1)
                        fprintf(['Qi and Ri found to guarantee asympt. stability '...
                                 'of the global system \n']);
                else
                        fprintf(['WARNING: no Qi and Ri found such that asympt. stability '...
                                 'of the global system is guaranteed \n']);
                end
             end 
             for i= obj.activeDGU % recompute neighbors passive gains
                 neighbors_i = sort([i; neighbors(obj.NetGraph, i)]); 
                 Q_Ni{i} = blkdiag(Qi{neighbors_i});
             end
        end
        
    
        %-----------------------TRANSITION PHASE-------------------------------%
        %- Transition Phase using Online reconfigurable terminal ingredients -%
        function [X, U, lenSim, xs, us, alpha, SolverTime] = transitionPhase(x0, paramBefore,...
                                                         paramAfter, Qi, Ri, target, ADMM, regulation)
            N = 5; %Horizon 
            % compute steady state
            if ADMM
                disp("Using ADMM !");
                [xs, us, alpha, SolverTime] = transition_compute_ss_admm(horzcat(x0{:}), N, paramBefore,...
                                    paramAfter, target);
            elseif ~ADMM
                disp("Not using ADMM !");
                [xs, us, alpha] = transition_compute_ss(horzcat(x0{:}), N, paramBefore,...
                                    paramAfter, target);
            else
                disp("Error: Choose ADMM to be true or false");
                
            end
            disp("Steady-State Xs = "); disp(xs);
            disp("Steady-State Us = "); disp(us);
            % MPC REGULATION to steady state
            X{1} = horzcat(x0{:});
            U{1} = 0; lenSim = 0;
            clear regulation2ss
            clear regulation2ss_admm
            n = 1;
            if regulation
                while any(abs(X{n}(1,paramBefore.activeDGU) - xs(1,paramBefore.activeDGU)) > 5e-2) || ...
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
        end
        
        %- Transition Phase using offline reconfigurable terminal ingredients -%
        %- Delta formulation required 
        function [X, U, lenSim, xs, us, SolverTime] = transitionPhaseDeltaADMM(x0, paramBefore,...
                                                         paramAfter, Qi, Ri, target, alpha, regulation)
            N = 5; %Horizon  
            X{1} = horzcat(x0{:});
            U{1} = 0; lenSim = 0;
            dX{1} = X{1}-horzcat(paramAfter.Xref{:});
            [dXs, dUs, SolverTime] = transition_compute_delta_ss_admm(dX{1}, N, paramBefore,...
                                        paramAfter, target, alpha);
            xs = dXs + horzcat(paramAfter.Xref{:}); % reference for system after PnP
            us = dUs + horzcat(paramAfter.Uref{:});
            disp("Steady-State Xs = "); disp(xs);
            disp("Steady-State Us = "); disp(us);
            clear regulation2ss_admm
            n = 1;
            if regulation
                while any(abs(X{n}(1,paramBefore.activeDGU) - xs(1,paramBefore.activeDGU)) > 1e-1) || ...
                      any(abs(X{n}(2,paramBefore.activeDGU) - xs(2,paramBefore.activeDGU)) > 5e-2)    
                        % control input is of size nu x M 

                        dU{n} = regulation2ss_admm(dX{n}, N, paramBefore, dXs, dUs, Qi, Ri); % get first control input
                        U{n} = dU{n} + horzcat(paramAfter.Uref{:});
                        if isnan(U{n})
                            error("Input to apply to controller is Nan at iteration %d",n);
                        end
                        for i=paramBefore.activeDGU % Only modify the activated DGU to reach s-s
                            neighbors_i = sort([i; neighbors(paramBefore.NetGraph, i)]); % get neighbors
                            % create neighbor state vector comprising the state of subsystem i
                            % and the state of it's neighbors (concatenated)
                            x_Ni = reshape(dX{n}(:,neighbors_i),[],1); 
                            dX{n+1}(:, i) = paramBefore.A_Ni{i}*x_Ni + paramBefore.Bi{i}*dU{n}(:,i);
                            X{n+1}(:,i) = dX{n+1}(:,i) + paramAfter.Xref{i};
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
        end
                                                     
        %% ------------------------- MPC CONTROLLERS--------------------------%%
        %----- TRACKING MPC WITH RECONFIGURABLE TERMINAL INGREDIENTS ----------%
        function [X,U, alpha, solverTimeTotal] = mpc_DGU_tracking(controller, x0, length_sim,param, Q_Ni,...
                                          Ri)
                                      
            X = cell(1,length_sim+1); % state at each timestep
            U = cell(1, length_sim); % control input at each timestep
            alpha= cell(1, length_sim);
            %X{1} = horzcat(x0{:}); % initial state columns are subsystem i
            % initialize all steps and states, otherwise yalmip returns NaN in
            % non-activated DGU, resulting in NaN in the states which are given
            % again to yalmip which does not support NaN as x0
            X(:) = {horzcat(x0{:})}; 
            controllerType = functions(controller);
            fprintf("--INFO: Using controller %s -- \n ", controllerType.function);
            clear (controllerType.file)
             
            N = 5; % Horizon
            solverTimeTotal = 0;
            for n = 1:length_sim % loop over all subsystem
                % control input is of size nu x M 
                fprintf("---- Simulation step %d ---- \n", n);
                [U{n}, alpha{n}, Tk] = controller(X{n}, Q_Ni, Ri, N, param); % get first control input
                solverTimeTotal = solverTimeTotal + Tk;
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
            fprintf("Total solver time for simulation time = %d \n", solverTimeTotal);
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
        function [X,U,alphaEvolution, solverTimeTotal] = mpc_sim_DGU_delta(controller, x0, length_sim, param,...
                                          alpha, Q_Ni, Ri, Gamma_Ni)
                                      
            dX = cell(1,length_sim+1); % state at each timestep
            dU = cell(1,length_sim); % control input at each timestep
            X = cell(1,length_sim+1); % state at each timestep
            U = cell(1,length_sim); % control input at each timestep
            alphaEvolution= cell(1, length_sim);
            % initialize all steps and states, otherwise yalmip returns NaN in
            % non-activated DGU, resulting in NaN in the states which are given
            % again to yalmip which does not support NaN as x0
            X(:) = {horzcat(x0{:})}; 
            dX(:) = {horzcat(x0{:})-horzcat(param.Xref{:})};

            N = 5; % Horizon
            controllerType = functions(controller);
            fprintf("--INFO: Using controller %s -- \n ", controllerType.function);
            clear (controllerType.file)
            solverTimeTotal = 0;
            for n = 1:length_sim % loop over all subsystem
                % control input is of size nu x M
               fprintf("---- Simulation step %d ---- \n", n);
               [dU{n}, dXend, Tk] = controller(dX{n}, alpha, Q_Ni, Ri, N, param); % get first control input
               solverTimeTotal = solverTimeTotal + Tk; 
               U{n} = dU{n} + horzcat(param.Uref{:});
                if isnan(dU{n})
                    error("Input to apply to controller is Nan at iteration %d",n);
                end
                for i=param.activeDGU
                    out_neighbors = [i; neighbors(param.NetGraph, i)]; % get neighbors
                    out_neighbors = sort(out_neighbors); % sorted neighbor list
                    % create neighbor state vector comprising the state of subsystem i
                    % and the state of it's neighbors (concatenated)
                    dx_Ni = reshape(dX{n}(:,out_neighbors),[],1); 
                    dXend_Ni = reshape(dXend(:,out_neighbors),[],1); 
                    % Apply first input to the system and get next state
                    dX{n+1}(:, i) = param.A_Ni{i}*dx_Ni + param.Bi{i}*dU{n}(:,i);
                    X{n+1}(:,i) = dX{n+1}(:,i) + param.Xref{i};
                    % update variable constraining terminal set of each subsystem
                    alpha(i) = alpha(i) + dXend_Ni'*Gamma_Ni{i}*dXend_Ni;
                end
                alphaEvolution{n} = alpha;
            end 
            fprintf("Total solver time for simulation time = %d \n", solverTimeTotal);
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
        
    end % end of static methods
end % class end