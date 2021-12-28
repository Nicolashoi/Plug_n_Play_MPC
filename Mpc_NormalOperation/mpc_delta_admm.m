function [u0, Xend,Tk] = mpc_delta_admm(x0, alpha, Q_Ni, Ri, N, param)
    rho = 0.25;
    Tk = 0; k = 2; l=1;
    centralStopCond = true;
    TMAX = 10;
    persistent localOptimizer
    if isempty(localOptimizer)
        for i = param.activeDGU
            localOptimizer{i} = init_optimizer(i, Q_Ni{i}, Ri{i}, N, param,rho);
        end
    end
   
    % Initialization   
    for i=param.activeDGU
        neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
        % Initialize the global copy z
        z_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*param.ni,N);
        
        % Initialize the Lagrange multipliers
        y_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*param.ni,N);
    end
    
    while(Tk < TMAX)
        elapsedTime = zeros(1,length(param.activeDGU));
        for i=param.activeDGU % for each subsystems
            % Solve a local optimizati% Loop while terminal terminal time not overrunnedon problem for each subsystem i
            [w_Ni{i,k}, vi{i,k}, elapsedTime(i)] = local_optim(localOptimizer{i}, x0,...
                                                 z_Ni{i,l}, y_Ni{i,l}, sqrt(alpha(i)));
      
            % Obtain the sorted list of neighbors of system i
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            for j = neighbors_i' % Update local copy of system i
                % estimation of neighbors j by system i
                wi{j,k,i}.xi = param.Wij{i}{j}*w_Ni{i,k}.x_Ni; 
            end
        end 
        Tk = Tk + max(elapsedTime);
        % Update global copy of each subsystem
        for i = param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            %each subsystem i averages it's state over the set of neighbors
            zi{i,k} = update_global_copy(wi(i,k,neighbors_i));
        end
        % Share the global copy with all the neighbors (update the Ni states)
        for i=param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            % Update the global vector zNi
            xi_cell =  cellfun(@(x) x.xi, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_Ni = vertcat(xi_cell{:});
            % Update the Lagrange Multipliers
            y_Ni_inter = diff_struct(w_Ni{i,k}, z_Ni{i,k});
            y_Ni{i,k} = add_struct(y_Ni{i,l}, ...
                            structfun(@(x) rho.*x, y_Ni_inter, 'Un', false)) ; 
        end
       
        % Terminal condition which is centralized (good for having an estimate 
        % of how much time iteration are needed
        if centralStopCond
            r_norm{k} = 0;
            s_norm{k} = 0;
            for i = param.activeDGU
                r_struct = diff_struct(w_Ni{i,k}, z_Ni{i,k});
                r_norm{k} = r_norm{k} + sum(vecnorm(r_struct.x_Ni,2));
                s_struct = diff_struct(z_Ni{i,k}, z_Ni{i,k-1});
                s_norm{k} = s_norm{k}+ N*rho^2*sum(vecnorm(s_struct.x_Ni,2));
            end
            if r_norm{k} < 0.1 && s_norm{k} < 0.1
                break;
            end
        end
        k = k+1;
        l = l+1;
    end
    fprintf(['Max elapsed time for a system to converge %d and '...
              'max number of iterations %d \n'], Tk, l-1);
    u0 = zeros(1,param.nb_subsystems);
    Xend = zeros(param.ni, param.nb_subsystems);
    for i=param.activeDGU
        u0(:,i) = vi{i,end}.ui(:,1);
        Xend(:,i) = wi{i,end,i}.xi(:,end);
    end
end

function [w_Ni, vi, elapsedTime] = local_optim(localOptimizer, x0,z_Ni,...
                                               y_Ni, alpha_i)
                        
    [solutionSet, ~, ~, ~, ~, optimTime] = localOptimizer(x0, z_Ni.x_Ni, y_Ni.x_Ni,...
                                                             alpha_i);
    w_Ni.x_Ni = solutionSet{2};
    vi.ui = solutionSet{1}; 
    elapsedTime= optimTime.solvertime;
end


function localOptimizer = init_optimizer(i, Q_Ni, Ri, N, param, rho)
    objective_i = 0;
    constraints_i = [];
    ni = param.ni;
    nu = param.nu;
    M = param.nb_subsystems;
    X0 = sdpvar(param.ni,M,'full'); % state as rows and system number as column
    n_Ni = size(param.A_Ni{i},2); % size of Neighbors set
    
    % variables as input to optimizer object
    z_Ni.x_Ni = sdpvar(n_Ni, N ,'full');   
    y_Ni.x_Ni = sdpvar(n_Ni, N ,'full');   
    alpha_i = sdpvar(1);
    % Constraint variables for the augmented Lagrangian
    eX_Ni_L = sdpvar(N,1,'full');  
    eX_Ni_Q = sdpvar(N,1,'full');  


    % Variables for Dynamics and constraints
    Xi = sdpvar(ni,N, 'full');
    Ui = sdpvar(nu,N-1, 'full');
    X_Ni = sdpvar(n_Ni, N, 'full');
    
    % For objective function
    Si = 1000*eye(param.ni);
    
    % obtain sorted list of neighbors of system i
    neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
    idx_Ni = logical(kron((neighbors_i==i), ones(param.ni,1)));
    
    % Initial condition
    constraints_i = [constraints_i, Xi(:,1) == X0(:,i)];
    constraints_i = [constraints_i, X_Ni(:,1) == reshape(X0(:,neighbors_i), [],1)];
   
    %% Planning Horizon Loop
    for n = 1:N-1 
        % Distributed Dynamics
        [constraints_i, objective_i] = dynamicsConstraints(constraints_i, objective_i,...
                                       n, i, param, idx_Ni, Xi, X_Ni,Ui, ....
                                       Q_Ni, Ri);
         % Constraints for augmented Lagrangian     
        constraints_i = [constraints_i, eX_Ni_L(n) >= y_Ni.x_Ni(:,n)'*...
                        (X_Ni(:,n)-z_Ni.x_Ni(:,n)), eX_Ni_Q(n) >= rho/2 * ...
                        (X_Ni(:,n)-z_Ni.x_Ni(:,n))'*(X_Ni(:,n)-z_Ni.x_Ni(:,n))];
    end
    constraints_i = [constraints_i, X_Ni(idx_Ni,N) == Xi(:,N)]; 
    % Horizon 1-> N constraints for Lagrangian, added to objective function
    constraints_i = [constraints_i, eX_Ni_L(N) >= y_Ni.x_Ni(:,N)'*...
                    (X_Ni(:,N)-z_Ni.x_Ni(:,N)), eX_Ni_Q(N) >= rho/2 * ...
                    (X_Ni(:,N)-z_Ni.x_Ni(:,N))'*(X_Ni(:,N)-z_Ni.x_Ni(:,N))];
    objective_i = objective_i + sum(eX_Ni_L) + sum(eX_Ni_Q);

    %% Terminal Set constraints
    % Add cost if we deviate from equilibrium state and if equilibrium state
    % deviates from the reference
     objective_i = objective_i + Xi(:,end)'*param.Pi{i}*Xi(:,end);
    %----- Terminal set condition (Reconfigurable Terminal Ingredients)--------% 
    LMI_terminal = [inv(param.Pi{i})*alpha_i, Xi(:,end);...
                        Xi(:,end)', alpha_i];
     constraints_i = [constraints_i, LMI_terminal >= 0];
    %--------------------------------------------------------------------------%
   
    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
    % Parameters to give to solver (updated global copy and lagrangians)
    parameters_in = {X0, z_Ni.x_Ni, y_Ni.x_Ni, alpha_i};
    solutions_out = {Ui, X_Ni};
    % Create optimizer object
    localOptimizer = optimizer(constraints_i,objective_i,ops,parameters_in,solutions_out);
   
end


function [constraints_i, objective_i] = dynamicsConstraints(constraints_i,objective_i,...
                                        n, i, param, idx_Ni,Xi,X_Ni,Ui,Q_Ni, Ri)

    constraints_i = [constraints_i, Xi(:,n+1) == param.A_Ni{i}*X_Ni(:,n)+...
                                                 param.Bi{i}*Ui(:,n)];
    constraints_i = [constraints_i, X_Ni(idx_Ni,n) == Xi(:,n)];                                       
    % State and input constraints
    constraints_i = [constraints_i, param.Gx_Ni{i} * X_Ni(:,n)...
                              <= param.fx_Ni{i}];
    constraints_i = [constraints_i, param.Gu_i{i} * Ui(:,n)...
                               <= param.fu_i{i}];
    objective_i = objective_i + X_Ni(:,n)'*Q_Ni*X_Ni(:,n) + Ui(:,n)'*Ri*Ui(:,n);
end


function zi = update_global_copy(wi)
    fn = fieldnames(wi{1});
    for ii = 1:numel(fn)
        % extract fieldname structure for every neighbors
        extract = cellfun(@(x) x.(fn{ii}), wi(:), 'Un', false); 
        % concatenate the same structure fieldname into 3rd dimension
        zi.(fn{ii}) = mean(cat(3, extract{:}),3); % mean along 3rd dimension
    end
end

function result = diff_struct(w_Ni, z_Ni)
    fn = fieldnames(w_Ni);
    for i = 1:numel(fn)
        result.(fn{i}) = w_Ni.(fn{i}) - z_Ni.(fn{i});
    end
end

function result = add_struct(struct1, struct2)
    fn = fieldnames(struct1);
    for i = 1:numel(fn)
        result.(fn{i}) = struct1.(fn{i}) + struct2.(fn{i});
    end
end