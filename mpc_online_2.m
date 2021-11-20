function u0 = mpc_online_2(x0, Ki, Q_Ni, Ri, Pi, N, param)
    persistent mpc_optimizer
    % initialize controller, if not done already
    if isempty(mpc_optimizer)
        mpc_optimizer = init_optimizer(Ki, Q_Ni, Ri, Pi, N, param);
    end
    [u0, ~, ~, ~, ~, feasibility]= mpc_optimizer(x0);
    disp(feasibility.infostr);
    %u0 = mpc_optimizer(x0, alpha);
end
 
function mpc_optimizer = init_optimizer(Ki, Q_Ni, Ri, Pi, N, param)
    %param = param_2_DGU;
    M = param.number_subsystem;
    %% create variables for optimizer
    nx = size(param.Ai{1},1);
    nu = size(param.Bi{1},2); % size 1
    % Input cell array of size N-1, each cell is an array of size nu*M
    U = sdpvar(repmat(nu,1,N-1), repmat(M,1,N-1),'full');
    Ue = sdpvar(nu, M,'full');
    % State cell array of size N (timestep), each cell is an array of size nx*M
    X = sdpvar(repmat(nx,1,N),repmat(M,1,N),'full'); % contains state of each subsystem i
    Xe = sdpvar(nx,M,'full');
    X_Ni = cell(M,N-1); % cell array for neighbor states
    X_eNi = cell(M,1); 
    X0 = sdpvar(nx,M,'full'); % state as rows and system number as column
    alpha = sdpvar(M,1,'full');
    ci = sdpvar(M,1, 'full');
    di = sdpvar(M,1, 'full');
    lambda = cell(M,1);%sdpvar(M,1, 'full'); % lambda_i
    objective = 0;
    constraints = [];
    S = cell(M,1);
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
    for i=1:M % loop over all subsystems
        S{i} = eye(size(param.Ai{i},1));
        m_Ni = size(param.A_Ni{i},2); % get size of set of Neighbors
        % obtain sorted list of neighbors of system i
        neighbors = sort([i;successors(param.graph, i)]);
        X_eNi{i} = sdpvar(m_Ni,1,'full'); % neighbor equilibrium state i
        lambda{i} = sdpvar(m_Ni,1, 'full');
        % constraint for X_eNi =  concat of neighbor x_ei
        constraints = [constraints, X_eNi{i} == ...
                                    reshape(Xe(:,neighbors),[],1)];
        %% Equilibrium constraints
        constraints = [constraints, Xe(:,i) == param.A_Ni{i}*X_eNi{i} + ...
                                                param.Bi{i}*Ue(:,i)];
        constraints = [constraints, Ue(:,i) == Ki{i}*X_eNi{i} + di(i)];  
       
        for k = 1:N-1 % Planning Horizon Loop
            %% Neighbor States for each ith sytem at the kth horizon iteration
            X_Ni{i,k} = sdpvar(m_Ni,1,'full'); % neighbor set of state i
            % add a constraint for the neighbor state i to be equal to the
            % concatenated subsystem neighbor state vectors
            constraints = [constraints, X_Ni{i,k} == ...
                                        reshape(X{k}(:,neighbors),[],1)];
           
            %% Distributed Dynamics
            constraints = [constraints, X{k+1}(:,i) == param.A_Ni{i}*X_Ni{i,k}+...
                                                       param.Bi{i}*U{k}(:,i)];
            %% State and input constraints
            constraints = [constraints, param.Gx_i{i} * X{k}(:,i)...
                                      <= param.fx_i{i}];
            constraints = [constraints, param.Gu_i{i} * U{k}(:,i)...
                                       <= param.fu_i{i}];
            %% Objective
            % sum of local functions l(xi,uf) %
            objective = objective + ...
                        (X_Ni{i,k}-X_eNi{i})'*Q_Ni{i}*(X_Ni{i,k}-X_eNi{i})+...
                        (U{k}(:,i)-Ue(:,i))'*Ri{i}*(U{k}(:,i)-Ue(:,i));  
                             
        end
        % Terminal cost
        objective = objective + ... 
                    (X{end}(:,i)-Xe(:,i))'*Pi{i}*(X{end}(:,i)-Xe(:,i))+...
                    (Xe(:,i) - param.Xref{i})'*S{i}*(Xe(:,i) - param.Xref{i});
        % Terminal Set condition
        constraints = [constraints, (X{end}(:,i)-ci(i))'*Pi{i}*(X{end}(:,i)-ci(i))...
                                    <= alpha(i)^2];                             
    end    
    % parameter for initial condition
    constraints = [constraints, X{1} == X0];
    
    %% Create optimizer object 
    ops = sdpsettings('verbose',1); %options
    parameters_in = {X0};
    %solutions_out = {[U{1}], [X{:}], [X_Ni{:}]}; % general form 
    solutions_out = U{1}; % get U0 for each subsystem, size nu x M
    mpc_optimizer = optimizer(constraints,objective,ops,parameters_in,solutions_out);
end