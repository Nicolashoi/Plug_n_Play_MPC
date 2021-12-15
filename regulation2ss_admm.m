function u0 = regulation2ss_admm(x0, N, param, xs, us, Qi, Ri)
    rho = 0.25;
    TMAX = 40;
    Tk = 0; k = 2; l=1;
    for i=param.activeDGU
        neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
        z_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*param.ni,N);
        y_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*param.ni,N);
    end
    while(Tk < TMAX)
        for i=param.activeDGU % loop over all subsystems
           
            [w_Ni{i,k}, vi{i,k}, elapsedTime] = local_optim(i, x0, xs(:,i), us(:,i),...
                             Qi{i}, Ri{i}, N, param, z_Ni{i,l}, y_Ni{i,l}, rho, k);
            Tk = Tk + elapsedTime;
            % obtain sorted list of neighbors of system i
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            %             idx_Ni = logical(kron((neighbors_i==i), ones(param.ni,1)));
            for j = neighbors_i'
                wi{j,k,i}.xi = param.Wij{i}{j}*w_Ni{i,k}.x_Ni; % estimation of neighbors j by system i
            end
        end  
        
        for i = param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            zi{i,k} = update_global_copy(wi(i,k,neighbors_i));
        end
        
        for i=param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            % Update the global vector zNi
            xi_cell =  cellfun(@(x) x.xi, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_Ni = vertcat(xi_cell{:});

            % Update Lagrange Multipliers
            y_Ni_inter = diff_struct(w_Ni{i,k}, z_Ni{i,k});
            y_Ni{i,k} = add_struct(y_Ni{i,l}, ...
                            structfun(@(x) rho.*x, y_Ni_inter, 'Un', false)) ; 
        end
        % Terminal condition which is centralized (good for having an estimate 
        % of how much time iteration are needed
        r_norm{k} = 0;
        s_norm{k} = 0;
        for i = param.activeDGU
            r_struct = diff_struct(w_Ni{i,k}, z_Ni{i,k});
            r_norm{k} = r_norm{k} + sum(vecnorm(r_struct.x_Ni,2));
            s_struct = diff_struct(z_Ni{i,k}, z_Ni{i,k-1});
            s_norm{k} = s_norm{k}+ N*rho^2*sum(vecnorm(s_struct.x_Ni,2));
        end
        if r_norm{k} < 0.3 && s_norm{k} < 0.3
            break;
        end
        k = k+1;
        l = l+1;
    end
    fprintf("Time elapsed to converge %d and number of iterations %d \n", Tk, l-1);
    u0 = zeros(1,param.nb_subsystems);
    for i=param.activeDGU
        u0(:,i) = vi{i,end}.ui(:,1);
    end
end

function [w_Ni, vi, elapsedTime] = local_optim(i, x0, xs_i, us_i, Qi, Ri, N, param, z_Ni, y_Ni, rho, k)
    persistent localOptimizer
    if k==2
        localOptimizer{i} = init_optimizer(xs_i, us_i, i, N,Qi,Ri, param, rho);
    end
    [solutionSet, ~, ~, ~, ~, optimTime] = localOptimizer{i}(x0, z_Ni.x_Ni, y_Ni.x_Ni);
    w_Ni.x_Ni = solutionSet{1};
    vi.ui = solutionSet{2};
    elapsedTime= optimTime.solvertime;
end


function localOptimizer = init_optimizer(xs_i, us_i, i, N, Qi, Ri, param, rho)
    objective_i = 0;
    constraints_i = [];
    ni = param.ni;
    nu = param.nu;
    M = param.nb_subsystems;
    n_Ni = size(param.A_Ni{i},2); % size of Neighbors set 
      
    % variables as input to optimizer object
    X0 = sdpvar(ni,M,'full'); % state as rows and system number as column
    z_Ni.x_Ni = sdpvar(n_Ni, N, 'full');
    y_Ni.x_Ni = sdpvar(n_Ni, N ,'full');

    % Variables for 1st optimization part (DGU PnP active but disconnected from
    % the rest of the network
    X_Ni = sdpvar(n_Ni, N, 'full');
    Xi = sdpvar(ni,N, 'full');
    Ui = sdpvar(nu,N-1, 'full');
      
    %% CONSTRAINTS DYNAMICS AND OBJECTIVE
    neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
    idx_Ni = logical(kron((neighbors_i==i), ones(ni,1)));
    
    constraints_i = [constraints_i, Xi(:,1) == X0(:,i)];
    constraints_i = [constraints_i, X_Ni(:,1) == reshape(X0(:,neighbors_i), [],1)];
    
 
    % Planning Horizon Loop
    for n = 1:N-1 
        % Distributed Dynamics
        [constraints_i, objective_i] = dynamicsConstraints(constraints_i,objective_i,...
                                        n, i, param, idx_Ni, Xi, X_Ni,Ui,...
                                       xs_i, us_i, Qi, Ri);
       % augmented Lagrangian         
       objective_i = objective_i + y_Ni.x_Ni(:,n)'*(X_Ni(:,n)-z_Ni.x_Ni(:,n)) + ...
                     rho/2 * (X_Ni(:,n)-z_Ni.x_Ni(:,n))'*(X_Ni(:,n)-z_Ni.x_Ni(:,n));
    end
    constraints_i = [constraints_i, X_Ni(idx_Ni,N) == Xi(:,N)];
    
    objective_i = objective_i + y_Ni.x_Ni(:,N)'*(X_Ni(:,N)-z_Ni.x_Ni(:,N)) + ...
                   rho/2 * (X_Ni(:,N)-z_Ni.x_Ni(:,N))'*(X_Ni(:,N)-z_Ni.x_Ni(:,N));
  
    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options

    parameters_in = {X0, z_Ni.x_Ni, y_Ni.x_Ni};

    solutions_out = {X_Ni, Ui};
    localOptimizer = optimizer(constraints_i,objective_i,ops,parameters_in,solutions_out);
    %solutionSet = mpc_optimizer(x0);

end

function [constraints_i, objective_i] = dynamicsConstraints(constraints_i,objective_i,...
                                        n, i, param, idx_Ni, Xi, X_Ni, Ui,...
                                        xs_i, us_i, Qi, Ri)

    constraints_i = [constraints_i, Xi(:,n+1) == param.A_Ni{i}*X_Ni(:,n)+...
                                               param.Bi{i}*Ui(:,n)];
    constraints_i = [constraints_i, X_Ni(idx_Ni,n) == Xi(:,n)];                                       
    % State and input constraints
    constraints_i = [constraints_i, param.Gx_i{i} * Xi(:,n)...
                              <= param.fx_i{i}];
    constraints_i = [constraints_i, param.Gu_i{i} * Ui(:,n)...
                               <= param.fu_i{i}];
    objective_i = objective_i + ...
                    10*(Xi(:,n)-xs_i)'*Qi*(Xi(:,n)-xs_i)+...
                    100*(Ui(:,n)-us_i)'*Ri*(Ui(:,n)-us_i);

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