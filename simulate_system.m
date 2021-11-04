%% Simulate system
% Author:
%   Nicolas Hoischen
% BRIEF: % Simulate time steps of a system, given a controller in input and the
         % initial state
% Following Algorithm 2 of "Distributed Synthesis and Control of Constrained
% Linear Systems"

function [X,U] = simulate_system(controller, x0, length_sim, simulation, param,...
                                 alpha, Q_Ni, Ri, Pi, Gamma_Ni)
    switch simulation
        case "MPC"
            [X,U] = mpc_simulation(controller, x0, length_sim, param,...
                                          alpha, Q_Ni, Ri, Pi, Gamma_Ni);
        case "Passivity" % passivity only
            M = param.number_subsystem;
            Kpass = cell(1,M); D = cell(1,M); Ppass = cell(1,M); 
            Gamma_pass = cell(1,M);
            for i= 1:M
            [Kpass{i}, D{i}, Ppass{i}, Gamma_pass{i}] = controller_passivity(...
                                             param.Ai{1}, param.Bi{1},...
                                             param.Ci{1}, param.Fi{1}, ...
                                             param.L_tilde, param.global_sysd.C);
            end
            [X,U] = passive_controller_sim(x0, length_sim, Kpass,param);
        otherwise
            error("Simulation not implemented");
       
    end
end


function [X,U] = passive_controller_sim(x0, length_sim, Kpass,param)
    Kblock = blkdiag(Kpass{:});
    Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
    Cd = param.global_sysd.C;
    Acl = Ad+Bd*Kblock;
%     sysd = ss(Acl,Bd,Cd,[], -1);
%     disp("Is global system passive ?"); disp(isPassive(sysd));
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = x0;
    for n=1:length_sim
        X{n+1} = Acl*X{n};
        U{n} = Kblock*X{n};
    end
end


function [X, U]= mpc_simulation(controller, x0, length_sim, param,...
                                          alpha, Q_Ni, Ri, Pi, Gamma_Ni)
  
        M = param.number_subsystem; % number of subsystems
        X = cell(length_sim+1,1); % state at each timestep
        U = cell(length_sim,1); % control input at each timestep
        X{1} = x0; % initial state
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
                                