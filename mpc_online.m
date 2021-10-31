%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:

function [X,U] = mpc_online(x0, alpha, Q_Ni, Ri, Pi, Gamma_Ni)
    param = param_coupled_oscillator;
    length_sim = 30;
    M = param.number_subsystem;
    % X = cell(size(param.global_sysd.A,1),1,length_sim+1);
    % U = cell(size(param.global_sysd.B,2),1,length_sim);
    X = cell(length_sim+1,1);
    U = cell(length_sim,1);
    X{1} = x0;
    for n = 1:length_sim
        u0 = get_input(x0, alpha, Q_Ni, Ri, Pi);
        U{n} = u0;
        for i=1:M
            neighbors = [i, successors(param.graph, i)];
            neighbors = sort(neighbors);
            %x_Ni = X{n}(:,i);
            x_Ni = [];
            for j = 1:length(neighbors) % PROBLEM FOR ORDER OF NEIGHBORS HAS TO BE SAME AS STATE SPACE
                x_Ni = [x_Ni; X{n}(:,neighbors(j))];   % create vector of neighbor states
            end 
            X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*u0(:,i);
            alpha(i) = alpha(i) + x_Ni'*Gamma_Ni{i}*x_Ni;
            x0(:,i) = X{n+1}(:, i);
        end
    end
end

function u0 = get_input(x0, alpha, Q_Ni, Ri, Pi)

    param = param_coupled_oscillator;
    N = 5;
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
    X0 = sdpvar(nx,M,'full'); % state as rows and system number as column
    alpha_var = sdpvar(M,1,'full');
    % CONSTRAINTS FOR EACH TIME STEP, GLOBAL CENTRALIZED SYSTEM
    for i=1:M % FOR EACH SUBSYSTEM
        x_Ni = cell(1,N);
        neighbors = [i,successors(param.graph, i)];
        %neighbors = sort(neighbors); % Sorting here is causing troubles
        for k = 1:N-1 % Planning Horizon
            
            for j = 1:length(neighbors)
                x_Ni{k} = [x_Ni{k}; X{k}(:,neighbors(j))];   % create vector of neighbor states
            end
            constraints = [constraints, X{k+1}(:,i) == param.A_Ni{i}*x_Ni{k}+...
                                                       param.Bi{i}*U{k}(:,i)];
     
            % state constraints
            constraints = [constraints, param.Gx_i{i} * X{k+1}(:,i)...
                                        <= param.fx_i{i}];
            % input constraints
            constraints = [constraints, param.Gu_i{i} * U{k}(:,i)...
                                        <= param.fu_i{i}];
            % sum of local functions l(xi,uf) %
            objective = objective + x_Ni{k}'*Q_Ni{i}*x_Ni{k} + U{k}(:,i)'*Ri{i}*...   
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
    solutions_out = {[U{1}], [X{:}]}; % general form 
    MPC_optimizer = optimizer(constraints,objective,ops,parameters_in,solutions_out);
    % get U0 for each subsystem
    solutions = MPC_optimizer({x0, alpha});
    u0 = solutions{1}; % size nu x M, columns are for each subsystem
    
end