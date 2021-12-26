%% Algorithm 2: Online Distributed MPC
% Author:
%   Nicolas Hoischen
% BRIEF: % Online mpc controller function. Executed at every subsystem. 
% Following Algorithm 2 of "Distributed Synthesis and Control of Constrained
% Linear Systems"

function [u0, Xend] = mpc_delta(x0, alpha, Q_Ni, Ri, N, param)
    persistent mpc_optimizer
    % initialize controller, if not done already
    if isempty(mpc_optimizer)
        mpc_optimizer = init_optimizer(Q_Ni, Ri, N, param);
    end
    [sol, ~, ~, ~, ~, feasibility]= mpc_optimizer(x0, sqrt(alpha));
    if isequal(feasibility.problem,1)
        disp(feasibility.infostr);
    end
    u0 = sol{1};
    Xend = sol{2};
end
 
function mpc_optimizer = init_optimizer(Q_Ni, Ri, N, param)
    M = param.nb_subsystems;
    %M = length(param.activeDGU);
    %% create variables for optimizer
    nx = size(param.Ai{1},1);
    nu = size(param.Bi{1},2); 
    % Input cell array of size N-1, each cell is an array of size nu*M
    U = sdpvar(repmat(nu,1,N-1), repmat(M,1,N-1),'full');
    % State cell array of size N (timestep), each cell is an array of size nx*M
    X = sdpvar(repmat(nx,1,N),repmat(M,1,N),'full');
    X_Ni = cell(M,N-1); % cell array for neighbor states
    X0 = sdpvar(nx,M,'full'); % state as rows and system number as column
    alpha_var = sdpvar(M,1,'full');
    objective = 0;
    constraints = [];
    
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
    for i= param.activeDGU%1:M % loop over all subsystems
        constraints = [constraints, X{1}(:,i) == X0(:,i)];
        n_Ni = size(param.A_Ni{i},2); % get size of set of Neighbors
        % obtain sorted list of neighbors of system i
        out_neighbors = sort([i; neighbors(param.NetGraph, i)]);
       
       
        for k = 1:N-1 % Planning Horizon Loop
            X_Ni{i,k} = sdpvar(n_Ni,1,'full'); % neighbor state i
            % add a constraint for the neighbor state i to be equal to the
            % concatenated subsystem neighbor state vectors
            constraints = [constraints, X_Ni{i,k} == ...
                                        reshape(X{k}(:,out_neighbors),[],1)];
            % Distributed Dynamics
            constraints = [constraints, X{k+1}(:,i) == param.A_Ni{i}*X_Ni{i,k}+...
                                                       param.Bi{i}*U{k}(:,i)];
            % state constraints
            constraints = [constraints, param.Gx_i{i} * X{k}(:,i)...
                                      <= param.fx_i{i}];
%             % input constraints
            constraints = [constraints, param.Gu_i{i} * U{k}(:,i)...
                                       <= param.fu_i{i}];
            % sum of local functions l(xi,uf) %
            objective = objective + X_Ni{i,k}'*Q_Ni{i}*X_Ni{i,k} + U{k}(:,i)'*Ri{i}*...   
                                    U{k}(:,i);
        end
        % Terminal cost
        objective = objective + X{end}(:,i)'*param.Pi{i}*X{end}(:,i);
%         constraints = [constraints, X{end}(:,i)'*param.Pi{i}*X{end}(:,i)...
%                                   <= alpha_var(i)];
 
        LMI_terminalSet = [inv(param.Pi{i})*alpha_var(i), X{end}(:,i);...
                        X{end}(:,i)', alpha_var(i)];
        constraints = [constraints, LMI_terminalSet >= 0];
    end    
    % parameter for initial condition
%     constraints = [constraints, X{1} == X0];
    
    %% Create optimizer object 
    ops = sdpsettings('verbose',1); %options
    parameters_in = {X0, alpha_var};
    %solutions_out = {[U{1}], [X{:}], [X_Ni{:}]}; % general form 
    solutions_out = {U{1}, [X{end}]}; % get U0 for each subsystem, size nu x M
    mpc_optimizer = optimizer(constraints,objective,ops,parameters_in,solutions_out);
end