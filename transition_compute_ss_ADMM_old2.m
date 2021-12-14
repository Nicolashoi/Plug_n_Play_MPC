
function [xs, us, alpha] = transition_compute_ss_ADMM_old2(x0, N, paramBefore, paramAfter, target)
    rho = 0.25;
    TMAX = 100;
    Tk = 0; k = 2; l=1;
    unionDGU = union(paramBefore.activeDGU,paramAfter.activeDGU);
    for i= paramBefore.activeDGU
        neighbors_i = sort([i;neighbors(paramBefore.NetGraph, i)]);
        z_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*paramBefore.ni,N);
        z_Ni{i,1}.x_eNi = zeros(length(neighbors_i)*paramBefore.ni,1);
        y_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*paramBefore.ni,N);
        y_Ni{i,1}.x_eNi = zeros(length(neighbors_i)*paramBefore.ni,1);
    end
    
    for i= paramAfter.activeDGU
        neighbors_i = sort([i;neighbors(paramAfter.NetGraph, i)]);
        z_Ni{i,1}.x_Ni_mod = zeros(length(neighbors_i)*paramAfter.ni,N);
        z_Ni{i,1}.alpha_Ni = zeros(length(neighbors_i)*paramAfter.ni,1);
        z_Ni{i,1}.c_Ni = zeros(length(neighbors_i)*paramAfter.ni,1); 
        y_Ni{i,1}.x_Ni_mod = zeros(length(neighbors_i)*paramAfter.ni,N);
        y_Ni{i,1}.alpha_Ni = zeros(length(neighbors_i)*paramAfter.ni,1);
        y_Ni{i,1}.c_Ni = zeros(length(neighbors_i)*paramAfter.ni,1); 
    end
    while(Tk < TMAX)
        tStart  =tic;
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
        for i = unionDGU % loop over all subsystems
            
          [w_Ni{i,k}, vi{i,k}] = local_optim(i,k, x0, N, paramBefore, paramAfter,...
                                 z_Ni{i,l}, y_Ni{i,l}, rho, target);
            neighbors_i = sort([i;union(neighbors(paramBefore.NetGraph, i), ...
                                neighbors(paramBefore.NetGraph, i))]);
            for j = neighbors_i'
                wi{j,k,i}.xi = horzcat(paramBefore.Wij{i}{j}*w_Ni{i,k}.x_Ni,...
                                       paramAfter.Wij{i}{j}*w_Ni{i,k}.x_Ni_mod);% estimation of neighbors j by system i
                                   
                wi{j,k,i}.xei = paramBefore.Wij{i}{j}*w_Ni{i,k}.x_eNi;
                extract_alpha_i = paramAfter.Wij{i}{j}*w_Ni{i,k}.alpha_Ni;
                wi{j,k,i}.alpha_i = extract_alpha_i(1); %array was alpha*dim(ni)
                wi{j,k,i}.ci = paramAfter.Wij{i}{j}*w_Ni{i,k}.c_Ni;
            end
        end  
        
        for i = unionDGU
            neighbors_i = sort([i;union(neighbors(paramBefore.NetGraph, i), ...
                                neighbors(paramBefore.NetGraph, i))]);
            zi{i,k} = update_global_copy(wi(i,k,neighbors_i));
        end
        
        for i=paramBefore.activeDGU
            neighbors_i = sort([i;neighbors(paramBefore.NetGraph, i)]);
            xi_cellBefore = cellfun(@(x) x.xi(:,1:N), zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_Ni = vertcat(xi_cellBefore{:});
            xei_cell =  cellfun(@(x) x.xei, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_eNi = vertcat(xei_cell{:});
        end
        
        for i=paramAfter.activeDGU
            neighbors_i = sort([i;neighbors(paramAfter.NetGraph, i)]);
            xi_cellAfter = cellfun(@(x) x.xi(:,N+1:end), zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_Ni_mod = vertcat(xi_cellAfter{:});
            alpha_i_cell = cellfun(@(x) x.alpha_i*ones(paramAfter.ni,1), ...
                           zi(neighbors_i,k), 'Un', false); % format is alpha_i*dim(ni)
            z_Ni{i,k}.alpha_Ni = vertcat(alpha_i_cell{:});
            ci_cell = cellfun(@(x) x.ci, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.c_Ni = vertcat(ci_cell{:});
            % Update Lagrange Multipliers
            y_Ni_inter = diff_struct(w_Ni{i,k}, z_Ni{i,k});
            y_Ni{i,k} = add_struct(y_Ni{i,l}, ...
                            structfun(@(x) rho.*x, y_Ni_inter, 'Un', false)) ; 
        end
        
        Tk = Tk + toc(tStart);
        % HERE PUT CONDITION FOR TERMINATION
        r_norm{k} = 0;
        s_norm{k} = 0;
        for i = unionDGU
            r_struct = diff_struct(w_Ni{i,k}, z_Ni{i,k});
            r_norm{k} = r_norm{k} + sum(vecnorm(r_struct.x_Ni,2));
            s_struct = diff_struct(z_Ni{i,k}, z_Ni{i,k-1});
            s_norm{k} = s_norm{k}+ N*rho^2*sum(vecnorm(s_struct.x_Ni,2));
        end
        if r_norm{k} < 1 && s_norm{k} < 1
            break;
        end
%         fprintf("Iteration %d,  Time elapsed %d \n", l, Tk);
        k = k+1;
        l = l+1;
       
    end
    xs = zeros(paramBefore.ni, paramBefore.nb_subsystems);
    us = zeros(paramBefore.nu, paramBefore.nb_subsystems);
    for i=unionDGU
        xs(:,i) = wi{i,end,i}.xei;
        us(:,i) = vi{i,end}.uei;
    end
end

function [w_Ni, vi] = local_optim(i,k, x0, N, paramBefore, paramAfter,z_Ni, y_Ni,...
                                  rho, target)
    persistent localOptimizer
    if k==2
        localOptimizer{i} = init_optimizer(x0,i, N,paramBefore, paramAfter,rho,...
                                            target);
    end
                                
    [solutionSet, ~, ~, ~, ~, elapsedtime] = localOptimizer{i}(z_Ni.x_Ni, z_Ni.x_eNi, z_Ni.alpha_Ni, z_Ni.c_Ni, ...
                     y_Ni.x_Ni, y_Ni.x_eNi, y_Ni.alpha_Ni, y_Ni.c_Ni);
    disp(elapsedtime.solvertime);
    w_Ni.x_Ni = solutionSet{4};
    w_Ni.x_Ni_mod = solutionSet{5};
    w_Ni.x_eNi = solutionSet{6};
    w_Ni.alpha_Ni = solutionSet{7};
    w_Ni.c_Ni = solutionSet{8};
    vi.ui = solutionSet{1};
    vi.uei = solutionSet{2};
    vi.di = solutionSet{3};

end

function localOptimizer = init_optimizer(x0,i, N, paramBefore, paramAfter, rho, target)
    objective_i = 0;
    constraints_i = [];
    ni = paramBefore.ni;
    nu = paramBefore.nu;
    M = paramBefore.nb_subsystems;
    Xi = sdpvar(ni,2*N, 'full');
    Ui = sdpvar(nu,2*N-1, 'full');
    n_Ni_before = size(paramBefore.A_Ni{i},2); % get size of set of Neighbors
    n_Ni_after = size(paramAfter.A_Ni{i},2); % get size of set of Neighbors
    X_eNi = sdpvar(n_Ni_before,1,'full'); % neighbor equilibrium state i
    X_Ni = sdpvar(n_Ni_before, N, 'full');
    Xei = sdpvar(ni,1,'full');
    Uei = sdpvar(nu,1,'full');
    di = sdpvar(nu,1,'full');
    
    
       % variables as input to optimizer object
    z_Ni.x_Ni = sdpvar(n_Ni_before, N, 'full');
    z_Ni.x_eNi = sdpvar(n_Ni_before, 1, 'full');
    z_Ni.x_Ni_mod = sdpvar(n_Ni_after, N, 'full');
    z_Ni.alpha_Ni = sdpvar(n_Ni_after,1, 'full');
    z_Ni.c_Ni = sdpvar(n_Ni_after,1, 'full');
    
    y_Ni.x_Ni = sdpvar(n_Ni_before, N ,'full');   
    y_Ni.x_eNi = sdpvar(n_Ni_before, 1, 'full');
    y_Ni.x_Ni_mod = sdpvar(n_Ni_after, N, 'full'); 
    y_Ni.alpha_Ni = sdpvar(n_Ni_after,1, 'full');
    y_Ni.c_Ni = sdpvar(n_Ni_after,1,'full'); 
    
    
    % Parameter for terminal set constraints
    X_Ni_mod = sdpvar(n_Ni_after, N, 'full');
    ci = sdpvar(ni,M,'full');
    alpha_i = sdpvar(1,M,'full');
    lambda_i = sdpvar(n_Ni_after,1,'full');
    bi = sdpvar(ni,1, 'full');
    %X0 = sdpvar(ni,M,'full'); % state as rows and system number as column
    constraints_i = [constraints_i, Xi(:,1) == x0(:,i)];
    % INCLUDE ??
    %constraints_i = [constraints_i, X_Ni(:,1) == vertcat(x0{neighbors_i})];
    neighbors_i = sort([i;neighbors(paramBefore.NetGraph, i)]);
    idx_Ni = logical(kron((neighbors_i==i), ones(ni,1)));
    %% Equilibrium constraints
    constraints_i = [constraints_i, Xei == paramBefore.A_Ni{i}*X_eNi + ...
                                            paramBefore.Bi{i}*Uei];
    constraints_i = [constraints_i, X_eNi(idx_Ni)==Xei];
    % With new redesigned local passive feedback gains
    constraints_i = [constraints_i, Uei == paramAfter.Ki{i}*Xei + di];  
    %constraints_i = [constraints_i,  paramAfter.Gx_Ni{i}*X_eNi <= paramAfter.fx_Ni{i}];
    %% Planning Horizon Loop
    
    for n = 1:N-1 
        % Distributed Dynamics
        [constraints_i, objective_i] = dynamicsConstraintsBeforePnP(constraints_i,objective_i,...
                                        n, i, paramBefore, idx_Ni,Xi,X_Ni,Ui,...
                                        x0(:,i), Xei, target);
       % augmented Lagrangian         
       objective_i = objective_i + y_Ni.x_Ni(:,n)'*(X_Ni(:,n)-z_Ni.x_Ni(:,n)) + ...
                     rho/2 * (X_Ni(:,n)-z_Ni.x_Ni(:,n))'*(X_Ni(:,n)-z_Ni.x_Ni(:,n));
    end
    % Terminal steady state condition for 1st optimization part
    constraints_i = [constraints_i, Xi(:,N) == Xei]; 
    constraints_i = [constraints_i, X_Ni(idx_Ni,N) == Xi(:,N)];
    % Augmented Lagrangian
    objective_i = objective_i + y_Ni.x_Ni(:,N)'*(X_Ni(:,N)-z_Ni.x_Ni(:,N)) + ...
                   y_Ni.x_eNi'*(X_eNi-z_Ni.x_eNi)+ ... 
                   rho/2 * (X_Ni(:,N)-z_Ni.x_Ni(:,N))'*(X_Ni(:,N)-z_Ni.x_Ni(:,N))...
                  + rho/2 * (X_eNi-z_Ni.x_eNi)'*(X_eNi-z_Ni.x_eNi);
    %% 2nd Optimization Part
    % Terminal Set constraints
    [constraints_i, alpha_Ni, c_Ni] = terminalConstraints(constraints_i, paramAfter,i,...
         ci, di,  alpha_i, lambda_i, bi);
    
    constraints_i = dynamicsConstraintsAfterPnP(constraints_i, ...
                                        objective_i, N, i, paramAfter, Xi, X_Ni_mod,...
                                        Ui, y_Ni, z_Ni, rho);
  
    % Terminal Set constraints
     constraints_i = [constraints_i, (Xi(:,end)-ci(:,i))'*paramAfter.Pi{i}*(Xi(:,end)-ci(:,i))...
                                   <= alpha_i(i)^2];
    % Augmented Lagrangian
    objective_i = objective_i + y_Ni.alpha_Ni'*(diag(alpha_Ni) - z_Ni.alpha_Ni) +...
                   y_Ni.c_Ni'*(c_Ni - z_Ni.c_Ni)+ rho/2*...
                  (diag(alpha_Ni) - z_Ni.alpha_Ni)'*(diag(alpha_Ni) - z_Ni.alpha_Ni)...
                  +rho/2 *(c_Ni - z_Ni.c_Ni)'*(c_Ni - z_Ni.c_Ni);
    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
   
    parameters_in = {z_Ni.x_Ni, z_Ni.x_eNi, z_Ni.alpha_Ni, z_Ni.c_Ni, ...
                     y_Ni.x_Ni, y_Ni.x_eNi,y_Ni.alpha_Ni, y_Ni.c_Ni};
    solutions_out = {Ui, Uei, di, X_Ni, X_Ni_mod, X_eNi, diag(alpha_Ni), c_Ni};
    localOptimizer = optimizer(constraints_i,objective_i,ops,parameters_in,solutions_out);
end


function [constraints_i, objective_i] = dynamicsConstraintsAfterPnP(constraints_i, ...
                                        objective_i, N, i, param, Xi, X_Ni_mod,...
                                        Ui, y_Ni, z_Ni, rho)
     
    % recompute neighbor set with parameters after Plug In / Plug Out
    neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
    idx_Ni = logical(kron((neighbors_i==i), ones(param.ni,1)));                                
    % X_Ni_mod is the neighbor state vector after Plug In /Plug out so shift indices for size of array   
     for n= N:2*N-1
        % Distributed Dynamics
        constraints_i = [constraints_i, X_Ni_mod(idx_Ni,n-N+1) == Xi(:,n)];  
        
        constraints_i = [constraints_i, Xi(:,n+1) == param.A_Ni{i}*X_Ni_mod(:,n-N+1)+...
                                                   param.Bi{i}*Ui(:,n)];
        constraints_i = [constraints_i, param.Gx_Ni{i} * X_Ni_mod(:,n-N+1)...
                              <= param.fx_Ni{i}];
       constraints_i = [constraints_i, param.Gu_i{i} * Ui(:,n)...
                                   <= param.fu_i{i}];   
                               
       objective_i = objective_i + y_Ni.x_Ni_mod(:,n-N+1)'*(X_Ni_mod(:,n-N+1)-...
                    z_Ni.x_Ni_mod(:,n-N+1))...
                     + rho/2 *(X_Ni_mod(:,n-N+1)-z_Ni.x_Ni_mod(:,n-N+1))'...
                         *(X_Ni_mod(:,n-N+1)-z_Ni.x_Ni_mod(:,n-N+1));
     end
                                                                        
                                    
end

function [constraints_i, objective_i] = dynamicsConstraintsBeforePnP(constraints_i,objective_i,...
                                        n, i, param, idx_Ni,Xi,X_Ni,Ui,...
                                        x0_i, Xei, target)

    constraints_i = [constraints_i, Xi(:,n+1) == param.A_Ni{i}*X_Ni(:,n)+...
                                               param.Bi{i}*Ui(:,n)];
    constraints_i = [constraints_i, X_Ni(idx_Ni,n) == Xi(:,n)];                                       
    % State and input constraints
    constraints_i = [constraints_i, param.Gx_i{i} * Xi(:,n)...
                              <= param.fx_i{i}];
    constraints_i = [constraints_i, param.Gu_i{i} * Ui(:,n)...
                               <= param.fu_i{i}];
    if target == "reference"
        objective_i = objective_i + 100*(Xei - param.Xref{i})'*(Xei - param.Xref{i});%1*norm(Xei-param.Xref{i},2);
    elseif target == "current state"
        objective_i = objective_i + norm(Xei-x0_i,2);
    else
        disp("objective not well defined, choose reference or current state");
    end
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
     constraints_i = [constraints_i, (param.A_Ni{i}+param.Bi{i}*param.K_Ni{i})...
                       *c_Ni + param.Bi{i}*di- ci(:,i) == bi];
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
                           param.Gu_i{i}(k,:)*di + sum_GxNorm2...
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