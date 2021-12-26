function u0 = regulation2ss(x0, N, param, xs, us, Qi, Ri)
    persistent mpc_optimizer
    % initialize controller, if not done already
    if isempty(mpc_optimizer)
        mpc_optimizer = init_optimizer(N, param, xs, us, Qi, Ri);
    end
    [u0, ~, ~, ~, ~, feasibility]= mpc_optimizer(x0);
    %disp(feasibility.infostr);
     %sol = mpc_optimizer(x0);
end
 
function mpc_optimizer = init_optimizer(N, param, xs, us, Qi, Ri)
    M = param.nb_subsystems;%length(param.activeDGU);
    %% create variables for optimizer
    nx = param.ni;
    nu = param.nu; % size 1
    % Input cell array of size N-1, each cell is an array of size nu*M
    U = sdpvar(repmat(nu,1,N-1), repmat(M,1,N-1),'full');
    % State cell array of size N (timestep), each cell is an array of size nx*M
    X = sdpvar(repmat(nx,1,N),repmat(M,1,N),'full'); % contains state of each subsystem i
    X0 = sdpvar(nx,M,'full'); % state as rows and system number as column
    %X0 = sdpvar(repmat(nx,1,M),ones(1,M), 'full');
    X_Ni = cell(M,N-1); % cell array for neighbor states of i
    objective = 0;
    constraints = [];
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
    for i=param.activeDGU % loop over all subsystems
        % loop over all active DGU but before connection or disconnection
        constraints = [constraints, X{1}(:,i) == X0(:,i)];
        ni = size(param.A_Ni{i},1); 
        n_Ni = size(param.A_Ni{i},2); % get size of set of Neighbors
        % obtain sorted list of neighbors of system i
        neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
        %% Planning Horizon Loop
        for n = 1:N-1 
            % Neighbor States for each ith sytem at the kth horizon iteration
            X_Ni{i,n} = sdpvar(n_Ni,1,'full'); % neighbor set of state i
            % add a constraint for the neighbor state i to be equal to the
            % concatenated subsystem neighbor state vectors
            constraints = [constraints, X_Ni{i,n} == ...
                                        reshape(X{n}(:,neighbors_i),[],1)];
           
            % Distributed Dynamics
            constraints = [constraints, X{n+1}(:,i) == param.A_Ni{i}*X_Ni{i,n}+...
                                                       param.Bi{i}*U{n}(:,i)];
            % State and input constraints
            constraints = [constraints, param.Gx_i{i} * X{n}(:,i)...
                                      <= param.fx_i{i}];
            constraints = [constraints, param.Gu_i{i} * U{n}(:,i)...
                                       <= param.fu_i{i}];
            % Objective
            objective = objective + ...
                        (X{n}(:,i)-xs(:,i))'*Qi{i}*(X{n}(:,i)-xs(:,i))+...
                        (U{n}(:,i)-us(:,i))'*Ri{i}*(U{n}(:,i)-us(:,i));                    
        end
        constraints = [constraints, X{end}(:,i) == xs(:,i)];
    end    
    % parameter for initial condition
    %constraints = [constraints, X{1} == X0];
    
    %% Create optimizer object 
    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
    parameters_in = {X0};
    solutions_out = U{1}; % general form 
    %solutions_out = {U{1}, U{2}, U{3}, U{4}}; % get U0 for each subsystem, size nu x M
    mpc_optimizer = optimizer(constraints,objective,ops,parameters_in,solutions_out);
end