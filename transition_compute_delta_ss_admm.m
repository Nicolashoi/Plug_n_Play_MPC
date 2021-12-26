
function [xs, us] = transition_compute_delta_ss_admm(x0, N, paramBefore, ...
                                                           paramAfter, target, alpha)
    rho = 0.25;
    TMAX = 100;
    Tk = 0; k = 2; l=1;
    % all the DGU concerned i.e. M U P (actual network + DGU's to be plugged in
    % or out
    unionDGU = union(paramBefore.activeDGU,paramAfter.activeDGU);
    % Initialization
    % DGUs active before PnP
    for i= paramBefore.activeDGU
        neighbors_i = sort([i;neighbors(paramBefore.NetGraph, i)]);
        z_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*paramBefore.ni,N);
        z_Ni{i,1}.x_eNi = zeros(length(neighbors_i)*paramBefore.ni,1);
        y_Ni{i,1}.x_Ni = zeros(length(neighbors_i)*paramBefore.ni,N);
        y_Ni{i,1}.x_eNi = zeros(length(neighbors_i)*paramBefore.ni,1);
    end
    % DGUs active after PnP
    for i= paramAfter.activeDGU
        neighbors_i = sort([i;neighbors(paramAfter.NetGraph, i)]);
        z_Ni{i,1}.x_Ni_mod = zeros(length(neighbors_i)*paramAfter.ni,N);
        % x_Ni_mod for states after plug in / out (horizon N+1 -> 2N)
        y_Ni{i,1}.x_Ni_mod = zeros(length(neighbors_i)*paramAfter.ni,N);
    end
    % Loop while terminal terminal time not overrunned
    while(Tk < TMAX)
        for i = unionDGU % loop over all subsystems
            [w_Ni{i,k}, vi{i,k}, elapsedTime] = local_optim(i,k, x0, N, paramBefore, paramAfter,...
                                 z_Ni{i,l}, y_Ni{i,l}, sqrt(alpha(i)), rho, target);
            Tk = Tk + elapsedTime; 
            % Obtain the sorted list of neighbors of system i: this list will be
            % used to update the local copies of each decision variable. The
            % intersection of the DGU neighbors connected before and after PnP
            % is taken, as we consider that an active DGU but disconnected from
            % the rest of the network can not estimate another state (as it is
            % not yet a neighbor/ or no longer a neighbor)
            
            neighbors_i = sort([i;intersect(neighbors(paramBefore.NetGraph, i), ...
                                neighbors(paramAfter.NetGraph, i))]);
                           
            for j = neighbors_i'% Update local copy of system i
                % Format wi(j,k,i) is estimation of system j by system i at
                % iteration k
                if any(paramAfter.activeDGU(:) == i) 
                    % Terminal ingredients for Horizon N->2N only for DGU active
                    % after PnP
                    wi{j,k,i}.xi = horzcat(paramBefore.Wij{i}{j}*w_Ni{i,k}.x_Ni,...
                                       paramAfter.Wij{i}{j}*w_Ni{i,k}.x_Ni_mod);
                else
                    % if DGU is not active after PnP (i.e. from horizon N -> 2N,
                    % we don't need to compute states after disconnection 
                    % set them by zero (default)
                    wi{j,k,i}.xi = horzcat(paramBefore.Wij{i}{j}*w_Ni{i,k}.x_Ni,...
                                   zeros(paramBefore.ni,N));
                end     
                wi{j,k,i}.xei = paramBefore.Wij{i}{j}*w_Ni{i,k}.x_eNi;
            end
        end  
        %% Update global copy of each subsystem
        for i = unionDGU
            % Again do the average estimate over the smallest subset (either
            % before PnP in which case a DGU is to be added (not yet connected),
            % or after PnP in which case a DGU is to be removed (disconnected).
            % As the average over Horizon 2N cannot be taken for disconnected
            % neighbors, the intersection (connected neighbors along the whole
            % procedure is taken
            neighbors_i = sort([i;intersect(neighbors(paramBefore.NetGraph, i), ...
                                neighbors(paramAfter.NetGraph, i))]);
            zi{i,k} = update_global_copy(wi(i,k,neighbors_i));
        end
        %% Share the global copy with all the neighbors (update the Ni states)
        for i=paramBefore.activeDGU
            neighbors_i = sort([i;neighbors(paramBefore.NetGraph, i)]);
            xi_cellBefore = cellfun(@(x) x.xi(:,1:N), zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_Ni = vertcat(xi_cellBefore{:});
            xei_cell =  cellfun(@(x) x.xei, zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_eNi = vertcat(xei_cell{:});
        end
        % variables at Horizon > N (for the new Network topology)
        for i=paramAfter.activeDGU
            neighbors_i = sort([i;neighbors(paramAfter.NetGraph, i)]);
            xi_cellAfter = cellfun(@(x) x.xi(:,N+1:end), zi(neighbors_i,k), 'Un', false);
            z_Ni{i,k}.x_Ni_mod = vertcat(xi_cellAfter{:});
        end
        for i=unionDGU
            % Update Lagrange Multipliers
            y_Ni_inter = diff_struct(w_Ni{i,k}, z_Ni{i,k});
            y_Ni{i,k} = add_struct(y_Ni{i,l}, ...
                            structfun(@(x) rho.*x, y_Ni_inter, 'Un', false)) ; 
        end
        %% Terminal condition which is centralized (good for having an estimate 
        % of how much time iteration are needed
        r_norm{k} = 0;
        s_norm{k} = 0;
        for i = unionDGU
            r_struct = diff_struct(w_Ni{i,k}, z_Ni{i,k});
            r_norm{k} = r_norm{k} + sum(vecnorm(r_struct.x_Ni,2));
            s_struct = diff_struct(z_Ni{i,k}, z_Ni{i,k-1});
            s_norm{k} = s_norm{k}+ N*rho^2*sum(vecnorm(s_struct.x_Ni,2));
        end
        if r_norm{k} < 0.025 && s_norm{k} < 0.025
            break;
        end
        fprintf("Iteration %d,  Time elapsed for each iteration %d \n", l, ...
                 elapsedTime);
        k = k+1;
        l = l+1;
       
    end
    %% Set variables
    xs = zeros(paramBefore.ni, paramBefore.nb_subsystems);
    us = zeros(paramBefore.nu, paramBefore.nb_subsystems);
    for i=unionDGU
        xs(:,i) = wi{i,end,i}.xei;
        us(:,i) = vi{i,end}.uei;
    end
    disp("Feasible steady-state found");
    fprintf("Total iterations %d and total time elapsed %d \n", l-1, Tk);
end

function [w_Ni, vi, elapsedTime] = local_optim(i,k, x0, N, paramBefore,...
                                   paramAfter,z_Ni, y_Ni, alpha_i, rho, target)
    persistent localOptimizer
    if k==2
        localOptimizer{i} = init_optimizer(x0,i, N,paramBefore, paramAfter,alpha_i,...
                                            rho, target);
    end
    
    if any(paramAfter.activeDGU(:) == i)
        [solutionSet, ~, ~, ~, ~, optimTime] = localOptimizer{i}(z_Ni.x_Ni, ...
         z_Ni.x_eNi, y_Ni.x_Ni, y_Ni.x_eNi);
                 
    w_Ni.x_Ni_mod = solutionSet{5};
    
    else
        [solutionSet, ~, ~, ~, ~, optimTime] = localOptimizer{i}(z_Ni.x_Ni, z_Ni.x_eNi, ...
                     y_Ni.x_Ni, y_Ni.x_eNi);
    end
                 
    if optimTime.problem
        error("Steady-state not found, Optimization problem: P&P rejected");
    end
    vi.ui = solutionSet{1};
    vi.uei = solutionSet{2};
    w_Ni.x_Ni = solutionSet{3};
    w_Ni.x_eNi = solutionSet{4};
     
    elapsedTime= optimTime.solvertime;
end

function localOptimizer = init_optimizer(x0,i, N, paramBefore, paramAfter,alpha_i,...
                                         rho, target)
    objective_i = 0;
    constraints_i = [];
    % Param
    ni = paramBefore.ni;
    nu = paramBefore.nu;
    M = paramBefore.nb_subsystems;
    % changing size parameters
    n_Ni_before = size(paramBefore.A_Ni{i},2); % size Neighbors's set former topology
    n_Ni_after = size(paramAfter.A_Ni{i},2); % size Neighbors's set new topology
    
    % variables as input to optimizer object
    % global copies
    z_Ni.x_Ni = sdpvar(n_Ni_before, N, 'full');
    z_Ni.x_eNi = sdpvar(n_Ni_before, 1, 'full');
    z_Ni.x_Ni_mod = sdpvar(n_Ni_after, N, 'full');
    
    % lagrangian
    y_Ni.x_Ni = sdpvar(n_Ni_before, N ,'full');   
    y_Ni.x_eNi = sdpvar(n_Ni_before, 1, 'full');
    y_Ni.x_Ni_mod = sdpvar(n_Ni_after, N, 'full'); 
     
    
       eX_Ni_L = sdpvar(N,1,'full');  
    eX_eNi_L = sdpvar(1, 1, 'full');
    eX_Ni_mod_L = sdpvar(N,1, 'full'); 
    eTerm_L = sdpvar(1,1, 'full');

    eX_Ni_Q = sdpvar(N,1,'full');  
    eX_eNi_Q = sdpvar(1, 1, 'full');
    eX_Ni_mod_Q = sdpvar(N,1, 'full'); 
    eTerm_Q = sdpvar(1,1, 'full');
    
    % States and Input of subsystem i
    Xi = sdpvar(ni,2*N, 'full');
    Ui = sdpvar(nu,2*N-1, 'full');
    % Variables for 1st optimization part, optimization to reach
    % equilibrium/steady state for the old network topology
    Uei = sdpvar(nu,1,'full');
    di = sdpvar(nu,1,'full');
    X_eNi = sdpvar(n_Ni_before,1,'full'); % neighbor equilibrium state i
    X_Ni = sdpvar(n_Ni_before, N, 'full');
    Xei = sdpvar(ni,1,'full');
    
    % Variables for 2nd optimization part, on new topology (system dynamics
    % influenced by neighbors
    X_Ni_mod = sdpvar(n_Ni_after, N, 'full');
    
    neighbors_i = sort([i;neighbors(paramBefore.NetGraph, i)]);
    idx_Ni = logical(kron((neighbors_i==i), ones(ni,1)));
    
    % Initial state constraints
    constraints_i = [constraints_i, Xi(:,1) == x0(:,i)];
    constraints_i = [constraints_i, X_Ni(:,1) == reshape(x0(:,neighbors_i), [],1)];
    
    %% Equilibrium constraints
    constraints_i = [constraints_i, Xei == paramBefore.A_Ni{i}*X_eNi + ...
                                            paramBefore.Bi{i}*Uei];
    constraints_i = [constraints_i, X_eNi(idx_Ni)==Xei]; 
    constraints_i = [constraints_i, eX_eNi_Q >= y_Ni.x_eNi'*(X_eNi-z_Ni.x_eNi), ...
                                    eX_eNi_L >= rho/2 * (X_eNi-z_Ni.x_eNi)'*...
                                    (X_eNi-z_Ni.x_eNi)];
    objective_i = objective_i + eX_eNi_L + eX_eNi_Q;
    %% Planning Horizon Loop 1->N for the 1st Optimization Part
    for n = 1:N-1 
        % Distributed Dynamics for old topology
        [constraints_i, objective_i] = dynamicsConstraintsBeforePnP(constraints_i,objective_i,...
                                        n, i, paramBefore, idx_Ni,Xi,X_Ni,Ui);
        constraints_i = [constraints_i, eX_Ni_L(n) >= y_Ni.x_Ni(:,n)'*(X_Ni(:,n)-z_Ni.x_Ni(:,n)),...
                                        eX_Ni_Q(n) >= rho/2 * (X_Ni(:,n)-z_Ni.x_Ni(:,n))'*(X_Ni(:,n)-z_Ni.x_Ni(:,n))];
       % augmented Lagrangian         
%        objective_i = objective_i + y_Ni.x_Ni(:,n)'*(X_Ni(:,n)-z_Ni.x_Ni(:,n)) + ...
%                      rho/2 * (X_Ni(:,n)-z_Ni.x_Ni(:,n))'*(X_Ni(:,n)-z_Ni.x_Ni(:,n));
    end
    
    % Terminal steady state condition for 1st optimization part
    constraints_i = [constraints_i, Xi(:,N) == Xei]; 
    constraints_i = [constraints_i, X_Ni(idx_Ni,N) == Xi(:,N)];
    
     % Augmented Lagrangian at Horizon N 
    constraints_i = [constraints_i, eX_Ni_L(N) >= y_Ni.x_Ni(:,N)'*(X_Ni(:,N)-z_Ni.x_Ni(:,N)),... 
                     eX_Ni_Q(N) >= rho/2 * (X_Ni(:,N)-z_Ni.x_Ni(:,N))'*(X_Ni(:,N)-z_Ni.x_Ni(:,N))];
    objective_i = objective_i + sum(eX_Ni_L) + sum(eX_Ni_Q);
%     Augmented Lagrangian
%     objective_i = objective_i + y_Ni.x_Ni(:,N)'*(X_Ni(:,N)-z_Ni.x_Ni(:,N)) + ...
%                    y_Ni.x_eNi'*(X_eNi-z_Ni.x_eNi)+ ... 
%                    rho/2 * (X_Ni(:,N)-z_Ni.x_Ni(:,N))'*(X_Ni(:,N)-z_Ni.x_Ni(:,N))...
%                   + rho/2 * (X_eNi-z_Ni.x_eNi)'*(X_eNi-z_Ni.x_eNi);
    
    if target == "reference"
        objective_i = objective_i + (Xei - 0)'*(Xei - 0);
    elseif target == "current state"
        objective_i = objective_i + (Xei-x0(:,i))'*(Xei-x0(:,i));
    else
        disp("objective not well defined, choose reference or current state");
    end
    
   
    %% 2nd Optimization Part N+1 -> 2N with new Network topology
    % Dynamics once steady state has been reached at horizon N
    if any(paramAfter.activeDGU(:) == i)
        constraints_i = dynamicsConstraintsAfterPnP(constraints_i, ...
                                        objective_i, N, i, paramAfter, Xi, X_Ni_mod,...
                                        Ui, y_Ni, z_Ni, eX_Ni_mod_L,eX_Ni_mod_Q,rho);
    
       % Terminal Set constraints
        LMI_terminal = [inv(paramAfter.Pi{i})*alpha_i, Xi(:,end);...
                        Xi(:,end)', alpha_i];
        constraints_i = [constraints_i, LMI_terminal >= 0];
%         constraints_i = [constraints_i, Xi(:,end)'*paramAfter.Pi{i}*Xi(:,end)...
%                                    <= alpha_i];          
    solutions_out = {Ui, Uei, X_Ni, X_eNi, X_Ni_mod};
    else
        solutions_out = {Ui, Uei, X_Ni, X_eNi};
    end
    parameters_in = {z_Ni.x_Ni, z_Ni.x_eNi,y_Ni.x_Ni, y_Ni.x_eNi};
    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
    localOptimizer = optimizer(constraints_i,objective_i,ops,parameters_in,solutions_out);
end

function [constraints_i, objective_i] = dynamicsConstraintsAfterPnP(constraints_i, ...
                                        objective_i, N, i, param, Xi, X_Ni_mod,...
                                        Ui, y_Ni, z_Ni, eX_Ni_mod_L,eX_Ni_mod_Q, rho)
     
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
                               
       constraints_i = [constraints_i, eX_Ni_mod_L(n-N+1) >= y_Ni.x_Ni_mod(:,n-N+1)'*(X_Ni_mod(:,n-N+1)-...
                        z_Ni.x_Ni_mod(:,n-N+1)), eX_Ni_mod_Q(n-N+1) >= rho/2 *...
                        (X_Ni_mod(:,n-N+1)-z_Ni.x_Ni_mod(:,n-N+1))'*...
                      (X_Ni_mod(:,n-N+1)-z_Ni.x_Ni_mod(:,n-N+1))];
                               
%        objective_i = objective_i + y_Ni.x_Ni_mod(:,n-N+1)'*(X_Ni_mod(:,n-N+1)-...
%                     z_Ni.x_Ni_mod(:,n-N+1))...
%                      + rho/2 *(X_Ni_mod(:,n-N+1)-z_Ni.x_Ni_mod(:,n-N+1))'...
%                          *(X_Ni_mod(:,n-N+1)-z_Ni.x_Ni_mod(:,n-N+1));
     end
     objective_i = objective_i + sum(eX_Ni_mod_L) + sum(eX_Ni_mod_Q);                                                                    
                                    
end

function [constraints_i, objective_i] = dynamicsConstraintsBeforePnP(constraints_i,objective_i,...
                                        n, i, param, idx_Ni,Xi,X_Ni,Ui)

    constraints_i = [constraints_i, Xi(:,n+1) == param.A_Ni{i}*X_Ni(:,n)+...
                                               param.Bi{i}*Ui(:,n)];
    constraints_i = [constraints_i, X_Ni(idx_Ni,n) == Xi(:,n)];                                       
    % State and input constraints
    constraints_i = [constraints_i, param.Gx_i{i} * Xi(:,n)...
                              <= param.fx_i{i}];
    constraints_i = [constraints_i, param.Gu_i{i} * Ui(:,n)...
                               <= param.fu_i{i}];
 
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