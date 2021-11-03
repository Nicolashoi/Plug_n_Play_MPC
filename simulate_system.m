%% Simulate system
% Author:
%   Nicolas Hoischen
% BRIEF: % Simulate time steps of a system, given a controller in input and the
         % initial state
% Following Algorithm 2 of "Distributed Synthesis and Control of Constrained
% Linear Systems"

function [X,U] = simulate_system(controller, x0, alpha, Q_Ni, Ri, Pi, Gamma_Ni, ...
                                 length_sim)
    param = param_coupled_oscillator;
    M = param.number_subsystem; % number of subsystems
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = x0; % initial state
    N = 10; % Horizon
    for n = 1:length_sim % loop over all subsystem
        % control input is of size nu x M
        U{n} = controller(X{n}, alpha, Q_Ni, Ri, Pi, N); % get first control input
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