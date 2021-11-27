%% Simulate system
% Author:
%   Nicolas Hoischen
% BRIEF: % Simulate time steps of a system, given a controller in input and the
         % initial state
% Following Algorithm 2 of "Distributed Synthesis and Control of Constrained
% Linear Systems"

function [X,U, Pinf] = simulate_system(controller, x0, length_sim, simulation, param,...
                                  Q, R, Pi, Gamma_Ni, alpha)
    Pinf = [];
    switch simulation
        case "MPC"
            [X,U] = mpc_sim_DGU_delta(controller, x0, length_sim, param,...
                                  alpha, Q, R, Pi, Gamma_Ni);
      case "MPC online"
        M = param.nb_subsystems;
        Kpass = cell(1,M); Ppass = cell(1,M); K_Ni = cell(1,M); 
        for i= 1:M
        [Kpass{i}, ~, Ppass{i}, ~] = controller_passivity(...
                                         param.Ai{i}, param.Bi{i},...
                                         param.Ci{i}, param.Fi{i}, ...
                                         param.L_tilde, param.global_sysd.C);
        disp("Passive controller Gain"); disp(Kpass{i});
        end
        for i=1:M
            neighbors = sort([i; successors(param.graph, i)]); 
            K_block = blkdiag(Kpass{neighbors}); %block matrix of all neighbors K
            K_Ni{i} = K_block(neighbors==i,:); % extract only row corresponding to subsystem i
        end
        [X,U] = mpc_DGU_tracking(controller, x0, length_sim, param,...
                                      K_Ni, Q, R, Ppass);   
                                  
        case "PASSIVITY" % passivity only
            M = param.nb_subsystems;
            Kpass = cell(1,M); D = cell(1,M); Ppass = cell(1,M); 
            Gamma_pass = cell(1,M);
            for i= 1:M
            [Kpass{i}, D{i}, Ppass{i}, Gamma_pass{i}] = controller_passivity(...
                                             param.Ai{i}, param.Bi{i},...
                                             param.Ci{i}, param.Fi{i}, ...
                                             param.L_tilde, param.global_sysd.C);
            disp("Passive controller Gain"); disp(Kpass{i});
            end
            [X,U] = sim_DGU_distributed(x0, length_sim, Kpass,param);
           
        case "LQR"
            Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
            [Klqr,Pinf, ~] = dlqr(Ad, Bd, Q,R);
            disp("LQR controller Gain"); disp(-Klqr);
            [X,U] = sim_global_DGU(x0, length_sim,-Klqr,param);

        otherwise
            error("Simulation not implemented");
       
    end
end

%% SIMULATION FOR DGU
function [X,U] = sim_DGU_distributed(x0, length_sim, K,param)
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
            dU{n}(:,i) = K{i}*dX{n}(:,i);
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


function [X,U] = sim_global_DGU(x0, length_sim, K,param)
    
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

function [X,U] = mpc_DGU_tracking(controller, x0, length_sim, param,...
                                          K, Q_Ni, Ri, Pi)
                                      
        M = param.nb_subsystems; % number of subsystems
        X = cell(length_sim+1,1); % state at each timestep
        U = cell(length_sim,1); % control input at each timestep
        X{1} = horzcat(x0{:}); % initial state columns are subsystem i    
   
        N = 10; % Horizon
        clear mpc_online_2
        for n = 1:length_sim % loop over all subsystem
            % control input is of size nu x M 
            U{n} = controller(X{n}, K, Q_Ni, Ri, Pi, N, param); % get first control input
            if isnan(U{n})
                error("Input to apply to controller is Nan at iteration %d",n);
            end
            for i=1:M
                neighbors = [i; successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(X{n}(:,neighbors),[],1); 
                X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
            end
        end  
end
