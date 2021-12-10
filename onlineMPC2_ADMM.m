function u0 = onlineMPC2_ADMM(x0, Ki, Q_Ni, Ri, Pi, N, param)
    persistent mpc_optimizer
    % initialize controller, if not done already
    if isempty(mpc_optimizer)
        mpc_optimizer = init_optimizer(Ki, Q_Ni, Ri, Pi, N, param);
    end
    [u0, ~, ~, ~, ~, feasibility]= mpc_optimizer(x0);
    %disp(feasibility.infostr);
     %sol = mpc_optimizer(x0);
end
 
function mpc_optimizer = init_optimizer(Ki, Q_Ni, Ri, Pi, N, param)
    M = param.nb_subsystems;%length(param.activeDGU);
    %% create variables for optimizer
    nx = param.ni;
    nu = param.nu; % size 1
    % Input cell array of size N-1, each cell is an array of size nu*M
    Ui = sdpvar(repmat(nu,1,N-1), repmat(M,1,N-1),'full');
    % State cell array of size N (timestep), each cell is an array of size nx*M
    Xi = sdpvar(repmat(nx,1,N),repmat(M,1,N),'full'); % contains state of each subsystem i
    X0 = sdpvar(nx,M,'full'); % state as rows and system number as column
    %X0 = sdpvar(repmat(nx,1,M),ones(1,M), 'full');
    X_Ni = cell(M,N-1); % cell array for neighbor states of i
    % Equilibrium state and input
    Ue = sdpvar(nu, M,'full');
    Xe = sdpvar(nx, M,'full');
    X_eNi = cell(M,1);  % equilibrium neighbor state
    % decision variables relative to constraints
    alpha = sdpvar(M,1,'full');
    ci = sdpvar(nx, M, 'full');
    di = sdpvar(nu, M, 'full');
    lambda = cell(M,1);
    bi = sdpvar(nx,M, 'full'); % for diagonal dominance
    % other variables
    c_Ni = cell(1,M);
    c_Ni(:) = {0}; % initialize all cells to zero
    alpha_Ni = cell(1,M);
    alpha_Ni(:) = {0};
    objective = 0;
  
    TMAX = 20;
    Tk = 0; k = 2; l=0;
     for i=param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            z_Ni{i,1}.x_Ni = zeros(length(neighbors_i),N);
            z_Ni{i,1}.x_eNi = zeros(length(neighbors_i),1);
            %z_Ni{i,1}.alpha_Ni = zeros(length(neighbors_i),1);
           % z_Ni{i,1}.c_Ni = 0; %%% TO CHECK SIZES

            y_Ni{i,1}.x_Ni = zeros(length(neighbors_i),N);
            y_Ni{i,1}.x_eNi = zeros(length(neighbors_i),1);
            %y_Ni{i,1}.alpha_Ni = zeros(length(neighbors_i),1);
            %y_Ni{i,1}.c_Ni = 0; %%% TO CHECK SIZES

    end
    while(Tk < TMAX)
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
        for i=param.activeDGU % loop over all subsystems
            Xi = sdpvar(param.ni,N, 'full');
            Ui = sdpvar(param.nu,N-1, 'full');
            n_Ni = size(param.A_Ni{i},2); % get size of set of Neighbors
            X_eNi = sdpvar(n_Ni,1,'full'); % neighbor equilibrium state i
            X_Ni = sdpvar(n_Ni, N, 'full');
            Xei = sdpvar(ni,1,'full');
            Uei = sdpvar(nu,1,'full');
            di = sdpvar(nu,1,'full');
            constraints_i = [constraints_i, Xi(:,1) == X0(:,i)];
            Si = 1000*eye(param.ni);

            % obtain sorted list of neighbors of system i
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);

            %lambda{i} = sdpvar(n_Ni,1, 'full');
            % constraint for X_eNi =  concat of neighbor x_ei
    %         constraints_i = [constraints_i, X_eNi == ...
    %                                     reshape(Xe(:,neighbors_i),[],1)];
            %% Equilibrium constraints
            constraints_i = [constraints_i, Xei == param.A_Ni{i}*X_eNi + ...
                                                    param.Bi{i}*Uei];
            constraints_i = [constraints_i, Uei == param.K_Ni{i}*X_eNi + di];  


            %% Planning Horizon Loop
            for n = 1:N-1 
                % Neighbor States for each ith sytem at the kth horizon iteration
    %             X_Ni{i,n} = sdpvar(n_Ni,1,'full'); % neighbor set of state i
                % add a constraint for the neighbor state i to be equal to the
                % concatenated subsystem neighbor state vectors
    %             constraints_i = [constraints_i, X_Ni{i,n} == ...
    %                                         reshape(X{n}(:,neighbors_i),[],1)];

                % Distributed Dynamics
                constraints_i = [constraints_i, Xi(:,n+1) == param.A_Ni{i}*X_Ni(:,n)+...
                                                           param.Bi{i}*Ui(:,n)];
                % State and input constraints
                constraints_i = [constraints_i, param.Gx_Ni{i} * X_Ni(:,n)...
                                          <= param.fx_Ni{i}];
    %             constraints = [constraints, param.Gx_i{i} * X{n}(:,i)...
    %                                       <= param.fx_i{i}];
                constraints_i = [constraints_i, param.Gu_i{i} * Ui(:,n)...
                                           <= param.fu_i{i}];
                % Objective
                objective = objective + ...
                            (X_Ni(:,n)-X_eNi)'*Q_Ni{i}*(X_Ni(:,n)-X_eNi)+...
                            (Ui(:,n)-Uei)'*Ri{i}*(Ui(:,n)-Uei) + ...
                            y_Ni{i,k}.x_Ni(:,n)*(w_Ni{i,k}.x_Ni(:,n)-z_Ni{i,k}.x_Ni{:,n)...
                            +;                    
            end
            %% Terminal cost
            objective = objective + y_Ni{i,k}.x_Ni*w_Ni{i,k}
%             objective = objective + ... 
%                         (Xi{end}(:,i)-Xe(:,i))'*Pi{i}*(Xi{end}(:,i)-Xe(:,i))+...
%                         (Xe(:,i) - param.Xref{i})'*Si*(Xe(:,i) - param.Xref{i});

            %%  Terminal Set condition
            constraints_i = [constraints_i, (Xi{end}(:,i)-ci(i))'*Pi{i}*(Xi{end}(:,i)-ci(i))...
                                        <= alpha(i)^2];      
            ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
            diagnostics = optimize(constraints, objective, ops);
            if diagnostics.problem == 1
               sprintf('MOSEK solver thinks it is infeasible for system %d', i)
            end
            wNi{i,k}.x_Ni = value(X_Ni);
            wNi{i,k}.x_eNi = value(X_eNi);
            wNi{i,k}.alpha_Ni = value(X_eNi);
            %wNi{i,k}.cNi = value(cNi);
            vi{i,k}.ui = value(Ui);
            vi{i,k}.uei = value(Uei);
            vi{i,k}.di = value(di);
            %vi{i,k}.lambda_ij = value(lambda_ij);
            wi{i,k,i}.xi = value(Xi);
            wi{i,k,i}.xei = value(Xei);
            %wi{i,k,i}.alpha_i = value(alpha_i); % ith value computed by system i
            %wi{i,k,i}.ci = value(ci);
            
        end    
        for i=param.activeDGU
            z{i,k}.xi = mean(horzcat(wi{i,k,:}.xi),2);
            z{i,k}.xei = mean(horzcat(wi{i,k,:}.xei),2); 
            %z{i,k}.alpha_i = mean(horzcat(wi{i,k,:}.alpha_i),2); 
            %z{i,k}.ci = mean(horzcat(wi{i,k,:}.ci),2); 
        end
        for i=param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            z_Ni{i,k}.x_Ni = vertcat(z{neighbors_i,k}.xi ,[],1);
            z_Ni{i,k}.x_eNi = vertcat(z{neighbors_i,k}.xei ,[],1);
            %z_Ni{i,k}.alpha_Ni = vertcat(z{neighbors_i,k}.alpha_i ,[],1);
            %z_Ni{i,k}.c_Ni = vertcat(z{neighbors_i,k}.ci ,[],1);
            
        end
    end
    % parameter for initial condition
    %constraints = [constraints, X{1} == X0];
    
    %% Create optimizer object 
    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
    parameters_in = {X0};
    %solutions_out = {[U{:}], [X_eNi{1}], [X_eNi{2}], [X_eNi{3}], di, Ue}; % general form 
    solutions_out = Ui{1}; % get U0 for each subsystem, size nu x M
    mpc_optimizer = optimizer(constraints_i,objective,ops,parameters_in,solutions_out);
end


function [w_Ni_new, vi_new] = lagrangian(N, p, constraints_i, objective_i, w_Ni, vi, z_Ni, y_Ni)
    fn = fieldname(w_Ni);
    for ii = 1:numel(fn)
        objective = objective + y_Ni.(fn{ii}).*(w_Ni.(fn{ii}) - z_Ni.(fn{ii})) + ...
                    p/2*(w_Ni.(fn{ii}) - z_Ni.(fn{ii})).^2;
    end
    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
       fprintf("MOSEK solver thinks it is infeasible");
    end
    for ii = 1:numel(fn)
        w_Ni_new.(fn{ii}) = value(w_Ni.(fn{ii}));
    end
    fnv = fieldname(vi);
    for iii = 1:numel(fnv)
        vi_new.(fnv{iii}) = value(vi.(fnv{iii}));
    end
end

