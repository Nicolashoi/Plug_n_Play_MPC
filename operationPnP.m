classdef operationPnP
    methods (Static)
        
        function obj = setPassiveControllers(obj)
            for i= obj.activeDGU
                [obj.Ki{i}, Di, obj.Pi{i}, Gamma_i] = controller_passivity(...
                                             obj.Ai{i}, obj.Bi{i},...
                                             obj.Ci{i}, obj.Fi{i}, ...
                                             obj.L_tilde, obj.global_sysd.C,i);
                disp("Passive controller Gain"); disp(obj.Ki{i});
                %disp(Gamma_i);
                disp(min(eig(Gamma_i)));
               
            end
            for i= obj.activeDGU
                neighbors_i = sort([i; neighbors(obj.NetGraph, i)]); 
                K_block = blkdiag(obj.Ki{neighbors_i}); %block matrix of all neighbors K
                obj.K_Ni{i} = K_block(neighbors_i==i,:); % extract only row corresponding to subsystem i
            end
        end
                
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
                disp(min(eig(Gamma_i)));
             end 
             for i= obj.activeDGU
                out_neighbors = sort([i; neighbors(obj.NetGraph, i)]); 
                K_block = blkdiag(obj.Ki{out_neighbors}); %block matrix of all neighbors K
                obj.K_Ni{i} = K_block(out_neighbors==i,:); % extract only row corresponding to subsystem i
             end
        end
        
    function [X,U] = sim_DGU_distributed(x0, length_sim, param)
        M = param.nb_subsystems; % number of subsystems
        X = cell(length_sim+1,1); % state at each timestep
        U = cell(length_sim,1); % control input at each timestep
        X{1} = horzcat(x0{:}); % columns are subsystem i

        if ~param.delta_config
            error("Switch to delta configuration to use Passivity to converge to ref");
        end

        dX = cell(length_sim+1,1); % state at each timestep
        dU = cell(length_sim,1); % control input at each timestep
        dX{1} = X{1} - horzcat(param.Xref{:});

        for n=1:length_sim
             for i=1:M
                dU{n}(:,i) = param.K{i}*dX{n}(:,i);
                U{n} = dU{n} + param.Uref{i};
                neighbors = [i; successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(dX{n}(:,neighbors),[],1); 
                % Apply first input to the system and get next state
                dX{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*dU{n}(:,i);
                X{n+1}(:,i) = dX{n+1}(:,i) + param.Xref{i};
             end
        end  
    end
    
    function [X, U, xs, us, alpha] = transitionPhase(x0, length_sim, paramBefore,...
                                                     paramAfter, Qi, Ri, target)
     [xs, us, alpha] = transition_compute_ss(horzcat(x0{:}), 10, paramBefore,...
                            paramAfter, target);
    
     X = cell(length_sim+1,1); % state at each timestep
     U = cell(length_sim,1); % control input at each timestep
        %X{1} = horzcat(x0{:}); % initial state columns are subsystem i
        % initialize all steps and states, otherwise yalmip returns NaN in
        % non-activated DGU, resulting in NaN in the states which are given
        % again to yalmip which does not support NaN as x0
     X(:) = {horzcat(x0{:})}; 
   
    N = 10; % Horizon
    clear regulation2ss
    for n = 1:length_sim % loop over all subsystem
            % control input is of size nu x M 
            U{n} = regulation2ss(X{n}, N, paramBefore, xs, us, Qi, Ri); % get first control input
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
    end  
         for i=1:paramBefore.nb_subsystems
            if ~any(i == paramBefore.activeDGU(:))
                X{1}(:,i) = [NaN;NaN];
                for n = 1:length_sim 
                    X{n+1}(:,i) = [NaN;NaN];
                    U{n}(:,i) = NaN;
                end
            end
        end                                   
                              
    end
    
    function [X,U] = mpc_sim_DGU_delta(controller, x0, length_sim, param,...
                                          alpha, Q_Ni, Ri, Pi, Gamma_Ni)
                                      
        M = param.nb_subsystems; % number of subsystems
        dX = cell(length_sim+1,1); % state at each timestep
        dU = cell(length_sim,1); % control input at each timestep
        X = cell(length_sim+1,1); % state at each timestep
        U = cell(length_sim,1); % control input at each timestep
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
                neighbors = [i; successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(dX{n}(:,neighbors),[],1); 
                % Apply first input to the system and get next state
                dX{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*dU{n}(:,i);
                X{n+1}(:,i) = dX{n+1}(:,i) + param.Xref{i};
                % update variable constraining terminal set of each subsystem
                alpha(i) = alpha(i) + x_Ni'*Gamma_Ni{i}*x_Ni;
            end
        end                                                                             
end

function [X,U] = mpc_DGU_tracking(controller, x0, length_sim,param, Q_Ni, Ri)
                                      
        X = cell(length_sim+1,1); % state at each timestep
        %X(:) = {[0;0]};
        U = cell(length_sim,1); % control input at each timestep
        %X{1} = horzcat(x0{:}); % initial state columns are subsystem i
        % initialize all steps and states, otherwise yalmip returns NaN in
        % non-activated DGU, resulting in NaN in the states which are given
        % again to yalmip which does not support NaN as x0
        X(:) = {horzcat(x0{:})}; 
   
        N = 10; % Horizon
        clear mpc_online_2
        for n = 1:length_sim % loop over all subsystem
            % control input is of size nu x M 
            U{n} = controller(X{n}, param.K_Ni, Q_Ni, Ri, param.Pi, N, param); % get first control input
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
            if ~any(i == param.activeDGU(:))
                X{1}(:,i) = [NaN;NaN];
                for n = 1:length_sim 
                    X{n+1}(:,i) = [NaN;NaN];
                    U{n}(:,i) = NaN;
                end
            end
        end
    end
    
    end
end