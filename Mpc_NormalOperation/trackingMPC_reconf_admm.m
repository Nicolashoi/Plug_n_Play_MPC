function [u0, alpha, Tk] = trackingMPC_reconf_admm(x0, Q_Ni, Ri, N, param)
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
        z_Ni{i,1}.x_eNi = zeros(length(neighbors_i)*param.ni,1);
        z_Ni{i,1}.alpha_Ni = zeros(length(neighbors_i)*param.ni,1);
        z_Ni{i,1}.c_Ni = zeros(length(neighbors_i)*param.ni,1); 
        % Initialize the Lagrange multipliers
        y_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*param.ni,N);
        y_Ni{i,1}.x_eNi = zeros(length(neighbors_i)*param.ni,1);
        y_Ni{i,1}.alpha_Ni = zeros(length(neighbors_i)*param.ni,1);
        y_Ni{i,1}.c_Ni = zeros(length(neighbors_i)*param.ni,1); 

    end
    
    while(Tk < TMAX)
        elapsedTime = zeros(1,length(param.activeDGU));
        for i=param.activeDGU % for each subsystems
            % Solve a local optimizati% Loop while terminal terminal time not overrunnedon problem for each subsystem i
            [w_Ni{i,k}, vi{i,k}, elapsedTime(i)] = local_optim(localOptimizer{i}, x0,...
                                                 z_Ni{i,l}, y_Ni{i,l});
            % Obtain the sorted list of neighbors of system i
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            for j = neighbors_i' % Update local copy of system i
                % estimation of neighbors j by system i
                wi{j,k,i}.xi = param.Wij{i}{j}*w_Ni{i,k}.x_Ni; 
                wi{j,k,i}.xei = param.Wij{i}{j}*w_Ni{i,k}.x_eNi;
                extract_alpha_i = param.Wij{i}{j}*w_Ni{i,k}.alpha_Ni;
                wi{j,k,i}.alpha_i = extract_alpha_i(1); %array was alpha*dim(ni)
                wi{j,k,i}.ci = param.Wij{i}{j}*w_Ni{i,k}.c_Ni;
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
            xei_cell =  cellfun(@(x) x.xei, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_eNi = vertcat(xei_cell{:});
            alpha_i_cell = cellfun(@(x) x.alpha_i*ones(param.ni,1), ...
                           zi(neighbors_i,k), 'Un', false); % format is alpha_i*dim(ni)
            z_Ni{i,k}.alpha_Ni = vertcat(alpha_i_cell{:});
            ci_cell = cellfun(@(x) x.ci, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.c_Ni = vertcat(ci_cell{:});
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
            if r_norm{k} < 0.05 && s_norm{k} < 0.05
                break;
            end
        end
        k = k+1;
        l = l+1;
    end
    fprintf(['Max elapsed time for a system to converge %d and '...
              'max number of iterations %d \n'], Tk, l-1);
    u0 = zeros(1,param.nb_subsystems);
    alpha = zeros(param.nb_subsystems,1);
    c = zeros(param.ni,param.nb_subsystems);
    for i=param.activeDGU
        u0(:,i) = vi{i,end}.ui(:,1);
        alpha(i) = wi{i,end,i}.alpha_i;
        c(:,i) = wi{i,end,i}.ci;
    end
end

function [w_Ni, vi, elapsedTime] = local_optim(localOptimizer, x0, z_Ni, y_Ni)
                   
    [solutionSet, ~, ~, ~, ~, optimTime] = localOptimizer(x0, z_Ni.x_Ni, z_Ni.x_eNi,...
                                            z_Ni.alpha_Ni, z_Ni.c_Ni, y_Ni.x_Ni,...
                                            y_Ni.x_eNi, y_Ni.alpha_Ni, y_Ni.c_Ni);
    w_Ni.x_Ni = solutionSet{4};
    w_Ni.x_eNi = solutionSet{5};
    w_Ni.alpha_Ni = solutionSet{6};
    w_Ni.c_Ni = solutionSet{7};
    vi.ui = solutionSet{1};
    vi.uei = solutionSet{2};
    vi.di = solutionSet{3};    
    elapsedTime= optimTime.solvertime;
end


function localOptimizer = init_optimizer(i, Q_Ni, Ri, N, param, rho)
    objective_i = 0;
    constraints_i = [];
    ni = param.ni;
    nu = param.nu;
    M = param.nb_subsystems;
    
    n_Ni = size(param.A_Ni{i},2); % size of Neighbors set
    X0 = sdpvar(param.ni,M,'full'); % state as rows and system number as column
    % variables as input to optimizer object
    z_Ni.x_Ni = sdpvar(n_Ni, N ,'full');   
    z_Ni.x_eNi = sdpvar(n_Ni, 1, 'full');
    z_Ni.alpha_Ni = sdpvar(n_Ni,1, 'full');
    z_Ni.c_Ni = sdpvar(n_Ni,1,'full'); 
    
    y_Ni.x_Ni = sdpvar(n_Ni, N ,'full');   
    y_Ni.x_eNi = sdpvar(n_Ni, 1, 'full');
    y_Ni.alpha_Ni = sdpvar(n_Ni,1, 'full');
    y_Ni.c_Ni = sdpvar(n_Ni,1,'full'); 
    
    % Constraint variables for the augmented Lagrangian
    eX_Ni_L = sdpvar(N,1,'full');  
    eX_eNi_L = sdpvar(1, 1, 'full');
    eTerm_L = sdpvar(1,1, 'full');

    eX_Ni_Q = sdpvar(N,1,'full');  
    eX_eNi_Q = sdpvar(1, 1, 'full');
    eTerm_Q = sdpvar(1,1, 'full');

    % Variables for Dynamics and constraints
    Xi = sdpvar(ni,N, 'full');
    Ui = sdpvar(nu,N-1, 'full');
    X_eNi = sdpvar(n_Ni,1,'full'); % neighbor equilibrium state i
    X_Ni = sdpvar(n_Ni, N, 'full');
    Xei = sdpvar(ni,1,'full');
    Uei = sdpvar(nu,1,'full');
    di = sdpvar(nu,1,'full');
    
    % For terminal set
    ci = sdpvar(ni,M,'full');
    alpha_i = sdpvar(1,M,'full');
    lambda_i = sdpvar(n_Ni,1,'full');
    bi = sdpvar(ni,1, 'full');
   
    % For objective function
    Si = 1000*eye(param.ni);
    
    % obtain sorted list of neighbors of system i
    neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
    idx_Ni = logical(kron((neighbors_i==i), ones(param.ni,1)));
    
    % Initial condition
    constraints_i = [constraints_i, Xi(:,1) == X0(:,i)];
    constraints_i = [constraints_i, X_Ni(:,1) == reshape(X0(:,neighbors_i), [],1)];
    %% Equilibrium constraints
    constraints_i = [constraints_i, Xei == param.A_Ni{i}*X_eNi + ...
                                            param.Bi{i}*Uei];
    constraints_i = [constraints_i, X_eNi(idx_Ni)==Xei];
    constraints_i = [constraints_i, Uei == param.K_Ni{i}*X_eNi + di];
    % Augmented Lagrangian constraints added to objective
    constraints_i = [constraints_i, eX_eNi_Q >= y_Ni.x_eNi'*(X_eNi-z_Ni.x_eNi), ...
                                    eX_eNi_L >= rho/2 * (X_eNi-z_Ni.x_eNi)'*...
                                    (X_eNi-z_Ni.x_eNi)];
    objective_i = objective_i + eX_eNi_L+ eX_eNi_Q ;

    %% Planning Horizon Loop
    for n = 1:N-1 
        % Distributed Dynamics
        [constraints_i, objective_i] = dynamicsConstraints(constraints_i, objective_i,...
                                       n, i, param, idx_Ni, Xi, X_Ni,Ui, ....
                                       X_eNi, Uei, Q_Ni, Ri);
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
    [constraints_i, alpha_Ni, c_Ni] = terminalConstraints(constraints_i, param,i,...
                                      ci, di,  alpha_i, lambda_i, bi);
    % Add cost if we deviate from equilibrium state and if equilibrium state
    % deviates from the reference
     objective_i = objective_i + ... 
                    (Xi(:,end)-Xei)'*param.Pi{i}*(Xi(:,end)-Xei)+...
                    (Xei - param.Xref{i})'*Si*(Xei - param.Xref{i});
    %----- Terminal set condition (Reconfigurable Terminal Ingredients)--------% 
%     LMI_terminal = [inv(param.Pi{i})*alpha_i(i), Xi(:,end) - ci(:,i);...
%                         Xi(:,end)' - ci(:,i)', alpha_i(i)];
%      constraints_i = [constraints_i, LMI_terminal >= 0];
     constraints_i = [constraints_i, norm(param.Pi{i}^(1/2)*(Xi(:,end)-ci(:,i)),2) <= alpha_i(i)];
    constraints_i = [constraints_i, alpha_i(i) >= 0];
    %--------------------------------------------------------------------------%
    % Constraints for augmented Lagrangian are added to objective
    constraints_i = [constraints_i, eTerm_L >= y_Ni.alpha_Ni'*(diag(alpha_Ni) ...
                    - z_Ni.alpha_Ni) + y_Ni.c_Ni'*(c_Ni - z_Ni.c_Ni),...
                    eTerm_Q >= rho/2*(diag(alpha_Ni) - z_Ni.alpha_Ni)'...
                    *(diag(alpha_Ni) - z_Ni.alpha_Ni) + rho/2 ...
                    *(c_Ni - z_Ni.c_Ni)'*(c_Ni - z_Ni.c_Ni)];
    objective_i = objective_i + eTerm_L +  eTerm_Q;  

    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
    % Parameters to give to solver (updated global copy and lagrangians)
    parameters_in = {X0, z_Ni.x_Ni, z_Ni.x_eNi, z_Ni.alpha_Ni, z_Ni.c_Ni, ...
                     y_Ni.x_Ni, y_Ni.x_eNi,y_Ni.alpha_Ni, y_Ni.c_Ni};
    solutions_out = {Ui, Uei, di, X_Ni, X_eNi, diag(alpha_Ni), c_Ni};
    % Create optimizer object
    localOptimizer = optimizer(constraints_i,objective_i,ops,parameters_in,solutions_out);
   
end


function [constraints_i, objective_i] = dynamicsConstraints(constraints_i,objective_i,...
                                        n, i, param, idx_Ni,Xi,X_Ni,Ui,...
                                        X_eNi, Uei, Q_Ni, Ri)

    constraints_i = [constraints_i, Xi(:,n+1) == param.A_Ni{i}*X_Ni(:,n)+...
                                                 param.Bi{i}*Ui(:,n)];
    constraints_i = [constraints_i, X_Ni(idx_Ni,n) == Xi(:,n)];                                       
    % State and input constraints
    constraints_i = [constraints_i, param.Gx_Ni{i} * X_Ni(:,n)...
                              <= param.fx_Ni{i}];
    constraints_i = [constraints_i, param.Gu_i{i} * Ui(:,n)...
                               <= param.fu_i{i}];
    objective_i = objective_i + ...
                    (X_Ni(:,n)-X_eNi)'*Q_Ni*(X_Ni(:,n)-X_eNi)+...
                    (Ui(:,n)-Uei)'*Ri*(Ui(:,n)-Uei);
end

function [constraints_i, alpha_Ni, c_Ni] = terminalConstraints(constraints_i, param,i,...
         ci, di,  alpha_i, lambda_i, bi)
     neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
     sumLambdaP_ij = 0;
     sumLambda_ij = 0;
     AbsSumLambdaP_ij = 0;
     alpha_Ni = 0;
     c_Ni = 0;
     
     for j = 1:length(neighbors_i)
         alpha_Ni = alpha_Ni + ...
                            alpha_i(neighbors_i(j))*(param.Wij{i}{neighbors_i(j)})'*... 
                            (param.Wij{i}{neighbors_i(j)});
         c_Ni = c_Ni + param.Wij{i}{neighbors_i(j)}'*ci(:,neighbors_i(j));
         Pij = (param.Wij{i}{neighbors_i(j)})'*param.Pi{neighbors_i(j)}...
                    *(param.Wij{i}{neighbors_i(j)});
         sumLambdaP_ij = sumLambdaP_ij + ...
                             lambda_i(j)*Pij;
         AbsSumLambdaP_ij = AbsSumLambdaP_ij + lambda_i(j)*abs(Pij);
         sumLambda_ij = sumLambda_ij + lambda_i(j);  
     end
     %% Equation 14: Approx of LMI with diagonal dominance
     PiInv = inv(param.Pi{i});
     constraints_i = [constraints_i, -bi <= (param.A_Ni{i}+param.Bi{i}*param.K_Ni{i})...
                       *c_Ni + param.Bi{i}*di- ci(:,i) <= bi];
     constraints_i = [constraints_i, alpha_i(i)-sumLambda_ij >= sum(bi(:))];  
     for k= 1:param.ni
            nondiag1 = sum(abs(PiInv(k,:))*alpha_i(i))- abs(PiInv(k,k))*alpha_i(i)+...
                       sum(abs(param.A_Ni{i}(k,:)+ param.Bi{i}(k)*param.K_Ni{i}(:))*alpha_Ni)...
                       + bi(k);
            constraints_i = [constraints_i, PiInv(k,k)*alpha_i(i) >= nondiag1];       
     end
     for k=1:size(param.A_Ni{i},2)
            nondiag2 = sum(AbsSumLambdaP_ij(k,:)) - AbsSumLambdaP_ij(k,k) + ...
                       sum(abs(param.A_Ni{i}(:,k)+ param.Bi{i}(:)*param.K_Ni{i}(k))*alpha_i(i));
            constraints_i = [constraints_i, sumLambdaP_ij(k,k) >= nondiag2];
     end
     
     %% Equation 11
        for k=1:size(param.Gx_Ni{i},1) 
            sum_GxNorm2 = 0;
            for j=1:length(neighbors_i)
                sum_GxNorm2 = sum_GxNorm2 + norm(param.Gx_Ni{i}(k,:)*...
                             (param.Wij{i}{neighbors_i(j)})'*...
                              param.Pi{neighbors_i(j)}^(-1/2),2) * alpha_i(neighbors_i(j));
                
            end
            constraints_i = [constraints_i, param.Gx_Ni{i}(k,:)*c_Ni + sum_GxNorm2 ...
                            <= param.fx_Ni{i}(k)];    
        end
        
        %% Equation 12
        for k=1:size(param.Gu_i{i},1) 
            sum_GuNorm2 = 0;
            for j=1:length(neighbors_i)
                sum_GuNorm2 = sum_GuNorm2 + norm(param.Gu_i{i}(k,:)*param.K_Ni{i}*...
                             (param.Wij{i}{neighbors_i(j)})'*...
                              param.Pi{neighbors_i(j)}^(-1/2),2) * alpha_i(neighbors_i(j));
                
            end
            constraints_i = [constraints_i, param.Gu_i{i}(k,:)*param.K_Ni{i}*c_Ni + ...
                           param.Gu_i{i}(k,:)*di + sum_GuNorm2...
                           <= param.fu_i{i}(k)];    
        end
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