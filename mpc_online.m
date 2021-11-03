%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:

function [X,U] = mpc_online(x0, alpha, Q_Ni, Ri, Pi, Gamma_Ni, length_sim)
    param = param_coupled_oscillator;
    M = param.number_subsystem; % number of subsystems
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = x0; % initial state
    
    for n = 1:length_sim % loop over all subsystem
        % control input is of size nu x M
        U{n} = get_input(X{n}, alpha, Q_Ni, Ri, Pi); % get first control input
        
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

function u0 = get_input(x0, alpha, Q_Ni, Ri, Pi)

    param = param_coupled_oscillator;
    N = 8;
    M = param.number_subsystem;
    % create variables for optimizer
    nx = size(param.Ai{1},1);
    nu = size(param.Bi{1},2); 
    objective = 0;
    constraints = [];
    % Input cell array of size N-1, each cell is an array of size nu*M
    U = sdpvar(repmat(nu,1,N-1), repmat(M,1,N-1),'full');
    % State cell array of size N (timestep), each cell is an array of size nx*M
    X = sdpvar(repmat(nx,1,N),repmat(M,1,N),'full');
    X_Ni = cell(M,N-1); % cell array for neighbor states
    X0 = sdpvar(nx,M,'full'); % state as rows and system number as column
    alpha_var = sdpvar(M,1,'full');
  
    for i=1:M % loop over all subsystems
        m_Ni = size(param.A_Ni{i},2); % get size of set of Neighbors
        % obtain sorted list of neighbors of system i
        neighbors = [i,successors(param.graph, i)];
        neighbors = sort(neighbors);
       
        for k = 1:N-1 % Planning Horizon Loop
            X_Ni{i,k} = sdpvar(m_Ni,1,'full'); % neighbor state i
            
            % add a constraint for the neighbor state i to be equal to the
            % concatenated subsystem neighbor state vectors
            constraints = [constraints, X_Ni{i,k} == ...
                                        reshape(X{k}(:,neighbors),[],1)];
            % Distributed Dynamics
            constraints = [constraints, X{k+1}(:,i) == param.A_Ni{i}*X_Ni{i,k}+...
                                                       param.Bi{i}*U{k}(:,i)];
     
            % state constraints
            constraints = [constraints, param.Gx_i{i} * X{k}(:,i)...
                                        <= param.fx_i{i}];
            % input constraints
            constraints = [constraints, param.Gu_i{i} * U{k}(:,i)...
                                        <= param.fu_i{i}];
            % sum of local functions l(xi,uf) %
            objective = objective + X_Ni{i,k}'*Q_Ni{i}*X_Ni{i,k} + U{k}(:,i)'*Ri{i}*...   
                                    U{k}(:,i);
        end
        % Terminal cost
        objective = objective + X{end}(:,i)'*Pi{i}*X{end}(:,i);
        constraints = [constraints, X{end}(:,i)'*Pi{i}*X{end}(:,i)...
                                    <= alpha_var(i)];
                                   
    end
              
    % parameter for initial condition
    constraints = [constraints, X{1} == X0];
    ops = sdpsettings('verbose',1); %options
    parameters_in = {X0, alpha_var};
    solutions_out = {[U{1}], [X{:}], [X_Ni{:}]}; % general form 
    MPC_optimizer = optimizer(constraints,objective,ops,parameters_in,solutions_out);
    % get U0 for each subsystem
    solutions = MPC_optimizer({x0, alpha});
    u0 = solutions{1}; % size nu x M, columns are for each subsystem
    
end