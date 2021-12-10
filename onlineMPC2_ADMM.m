function u0 = onlineMPC2_ADMM(x0,Q_Ni, Ri, N, param)
    p = 1/2;
    TMAX = 10;
    Tk = 0; k = 2; l=1;
     for i=param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            z_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*param.ni,N);
            z_Ni{i,1}.x_eNi = zeros(length(neighbors_i)*param.ni,1);
            %z_Ni{i,1}.alpha_Ni = zeros(length(neighbors_i),1);
           % z_Ni{i,1}.c_Ni = 0; %%% TO CHECK SIZES

            y_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*param.ni,N);
            y_Ni{i,1}.x_eNi = zeros(length(neighbors_i)*param.ni,1);
            %y_Ni{i,1}.alpha_Ni = zeros(length(neighbors_i),1);
            %y_Ni{i,1}.c_Ni = 0; %%% TO CHECK SIZES

    end
    while(Tk < TMAX)
        tStart  =tic;
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
        for i=param.activeDGU % loop over all subsystems
            
          [w_Ni{i,k}, vi{i,k}] = local_optim(i, x0, Q_Ni, Ri, N, param, z_Ni{i,l}, y_Ni{i,l});
            % obtain sorted list of neighbors of system i
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);

            for j = neighbors_i'
                wi{j,k,i}.xi = param.Wij{i}{j}*w_Ni{i,k}.x_Ni(:,:); % estimation of neighbors j by system i
                wi{j,k,i}.xei = param.Wij{i}{j}*w_Ni{i,k}.x_eNi;
            end
        end  
        
        for i = param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            zi{i,k} = update_global_copy(wi(i,k,neighbors_i));
        end
        for i=param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            xi_cell =  cellfun(@(x) x.xi, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_Ni = vertcat(xi_cell{:});
            xei_cell =  cellfun(@(x) x.xei, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_eNi = vertcat(xei_cell{:});
            y_Ni_inter = diff_struct(w_Ni{i,k}, z_Ni{i,k});
            y_Ni{i,k} = add_struct(y_Ni{i,l}, ...
                            structfun(@(x) p.*x, y_Ni_inter, 'Un', false)) ; 
        end
        Tk = Tk + toc(tStart);
        k = k+1;
        l = l+1;
    end
    u0 = zeros(1,param.nb_subsystems);
    for i=param.activeDGU
        u0(:,i) = vi{i,end}.ui(:,1);
    end
end


function [w_Ni, vi] = local_optim(i, x0, Q_Ni, Ri, N, param, z_Ni, y_Ni)
    neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
    p = 1/2;
    objective_i = 0;
    constraints_i = [];
    Xi = sdpvar(param.ni,N, 'full');
    Ui = sdpvar(param.nu,N-1, 'full');
    n_Ni = size(param.A_Ni{i},2); % get size of set of Neighbors
    X_eNi = sdpvar(n_Ni,1,'full'); % neighbor equilibrium state i
    X_Ni = sdpvar(n_Ni, N, 'full');
    Xei = sdpvar(param.ni,1,'full');
    Uei = sdpvar(param.nu,1,'full');
    di = sdpvar(param.nu,1,'full');
    ci = sdpvar(param.ni,1,'full');
    alpha = sdpvar(1);
    lambda = sdpvar(n_Ni,1,'full');
    bi = sdpvar(param.ni,1, 'full');
    
    constraints_i = [constraints_i, Xi(:,1) == x0{i}];
    for j=neighbors_i'
       constraints_i = [constraints_i, X_Ni(:,1) == vertcat(x0{neighbors_i})];
    end
    Si = 1000*eye(param.ni);

    % obtain sorted list of neighbors of system i
    neighbors_i = sort([i;neighbors(param.NetGraph, i)]);

    %% Equilibrium constraints
    constraints_i = [constraints_i, Xei == param.A_Ni{i}*X_eNi + ...
                                            param.Bi{i}*Uei];
    constraints_i = [constraints_i, Uei == param.K_Ni{i}*X_eNi + di];  

    
    %% Planning Horizon Loop
    for n = 1:N-1 
        % Distributed Dynamics
        constraints_i = [constraints_i, Xi(:,n+1) == param.A_Ni{i}*X_Ni(:,n)+...
                                                   param.Bi{i}*Ui(:,n)];
        idx_Ni = logical(kron((neighbors_i==i), ones(param.ni,1)));
        constraints_i = [constraints_i, X_Ni(idx_Ni,n+1) == Xi(:,n+1)];                                       
        % State and input constraints
        constraints_i = [constraints_i, param.Gx_Ni{i} * X_Ni(:,n)...
                                  <= param.fx_Ni{i}];
        constraints_i = [constraints_i, param.Gu_i{i} * Ui(:,n)...
                                   <= param.fu_i{i}];
        % Objective
        objective_i = objective_i + ...
                    (X_Ni(:,n)-X_eNi)'*Q_Ni{i}*(X_Ni(:,n)-X_eNi)+...
                    (Ui(:,n)-Uei)'*Ri{i}*(Ui(:,n)-Uei);
       % augmented Lagrangian         
       objective_i = objective_i + y_Ni.x_Ni(:,n)'*(X_Ni(:,n)-z_Ni.x_Ni(:,n)) + ...
                     p/2 * (X_Ni(:,n)-z_Ni.x_Ni(:,n))'*(X_Ni(:,n)-z_Ni.x_Ni(:,n));
    end
    objective_i = objective_i + y_Ni.x_Ni(:,end)'*(X_Ni(:,end)-z_Ni.x_Ni(:,end)) + ...
                   y_Ni.x_eNi'*(X_eNi-z_Ni.x_eNi)+ ... 
                   p/2 * (X_Ni(:,end)-z_Ni.x_Ni(:,end))'*(X_Ni(:,end)-z_Ni.x_Ni(:,end))...
                  + p/2 * (X_eNi-z_Ni.x_eNi)'*(X_eNi-z_Ni.x_eNi);
    %% Terminal cost
    objective_i = objective_i + ...
                (Xi(:,end)-Xei)'*param.Pi{i}*(Xi(:,end)-Xei) +...
                (Xei-param.Xref{i})'*Si*(Xei-param.Xref{i});
    
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0); %options
    diagnostics = optimize(constraints_i, objective_i, ops);
    if diagnostics.problem == 1
       fprintf("MOSEK solver thinks it is infeasible");
    end
    
    w_Ni.x_Ni = value(X_Ni);
    w_Ni.x_eNi = value(X_eNi);
    vi.ui = value(Ui);
    vi.uei = value(Uei);
    vi.di = value(di);
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