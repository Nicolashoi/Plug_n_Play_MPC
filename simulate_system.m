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
            [X,U] = mpc_simulation(controller, x0, length_sim, param,...
                                          alpha, Q, R, Pi, Gamma_Ni);
        case "PASSIVITY" % passivity only
            M = param.number_subsystem;
            Kpass = cell(1,M); D = cell(1,M); Ppass = cell(1,M); 
            Gamma_pass = cell(1,M);
            for i= 1:M
            [Kpass{i}, D{i}, Ppass{i}, Gamma_pass{i}] = controller_passivity(...
                                             param.Ai{1}, param.Bi{1},...
                                             param.Ci{1}, param.Fi{1}, ...
                                             param.L_tilde, param.global_sysd.C);
            end
            %disp("Passive controller Gain"); disp(Kpass{:});
            if param.name == "COUPLED_OSCI"
                [X,U] = sim_oscill_distributed(x0, length_sim, Kpass,param);
            elseif param.name == "2_DGU"
                [X,U] = sim_DGU_distributed(x0, length_sim, Kpass,param);
            end
        case "LQR"
            Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
            [Klqr,Pinf, ~] = dlqr(Ad, Bd, Q,R);
            disp("LQR controller Gain"); disp(-Klqr);
            if param.name == "COUPLED_OSCI"
                [X,U] = sim_global_coupled_osc(x0, length_sim, -Klqr,param);
            elseif param.name == "2_DGU"
                [X,U] = sim_global_DGU(x0, length_sim,-Klqr,param);
            end
        otherwise
            error("Simulation not implemented");
       
    end
end

function [X,U] = sim_oscill_distributed(x0, length_sim, K,param)
    M = param.number_subsystem; % number of subsystems
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = horzcat(x0{:}); % columns are subsystem i
    for n=1:length_sim
         for i=1:M
            U{n}(:,i) = K{i}*X{n}(:,i);
            neighbors = [i, successors(param.graph, i)]; % get neighbors
            neighbors = sort(neighbors); % sorted neighbor list
            % create neighbor state vector comprising the state of subsystem i
            % and the state of it's neighbors (concatenated)
            x_Ni = reshape(X{n}(:,neighbors),[],1); 
            % Apply first input to the system and get next state
            X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
         end
    end  
end

function [X,U] = sim_global_coupled_osc(x0, length_sim, K,param)
    
    Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
    %Bd(1:2,2) = [0;0]; Bd(3:4,1) = [0;0];
    Acl = Ad+Bd*K;
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = vertcat(x0{:}); % global initial state, vector Mx1
    for n=1:length_sim
        X{n+1} = Acl*X{n};
        U{n} = K*X{n};
    end     
end

function [X,U] = sim_DGU_distributed(x0, length_sim, K,param)
    M = param.number_subsystem; % number of subsystems
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = horzcat(x0{:}); % columns are subsystem i
    for n=1:length_sim
         for i=1:M
            X{n}(:,i) = X{n}(:,i) + [0;0;-param.Vr{i}];
            U{n}(:,i) = K{i}*X{n}(:,i) - param.Vr{i}/param.Vin{i} + ...
                        K{i}*[-param.Vr{i}; 0; 0];
            neighbors = [i, successors(param.graph, i)]; % get neighbors
            neighbors = sort(neighbors); % sorted neighbor list
            % create neighbor state vector comprising the state of subsystem i
            % and the state of it's neighbors (concatenated)
            x_Ni = reshape(X{n}(:,neighbors),[],1); 
            % Apply first input to the system and get next state
            X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
         end
    end  
end


function [X,U] = sim_global_DGU(x0, length_sim, K,param)
    
    Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
    M = param.number_subsystem; % number of subsystems
    Acl = Ad+Bd*K;
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = vertcat(x0{:}); % global initial state, vector Mx1
    for n=1:length_sim
        xconst = cell(1,M); uconst = cell(1,M); ufwd = cell(1,M);
        for i =1:M
            xconst{i} = [0;0;-param.Vr{i}];
            ufwd{i} = [-param.Vr{i}; 0; 0];
            uconst{i} = - param.Vr{i}/param.Vin{i}; 
        end
        X{n} = X{n} + vertcat(xconst{:});
        U{n} = K*X{n} + vertcat(uconst{:}) + K*vertcat(ufwd{:});
        X{n+1} = Ad*X{n} + Bd*U{n};
    end     
end


function [X, U]= mpc_simulation(controller, x0, length_sim, param,...
                                          alpha, Q_Ni, Ri, Pi, Gamma_Ni)
  
        M = param.number_subsystem; % number of subsystems
        X = cell(length_sim+1,1); % state at each timestep
        U = cell(length_sim,1); % control input at each timestep
        X{1} = horzcat(x0{:}); % initial state columns are subsystem i
        N = 10; % Horizon
        for n = 1:length_sim % loop over all subsystem
            % control input is of size nu x M
            U{n} = controller(X{n}, alpha, Q_Ni, Ri, Pi, N); % get first control input
            if isnan(U{n})
                error("Input to apply to controller is Nan at iteration %d",n);
            end
            for i=1:M
                neighbors = [i, successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(X{n}(:,neighbors),[],1); 
                % Apply first input to the system and get next state
                X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
                % update variable constraining terminal set of each subsystem
                alpha(i) = alpha(i) + x_Ni'*Gamma_Ni{i}*x_Ni;
            end
        end                                     
end
                                