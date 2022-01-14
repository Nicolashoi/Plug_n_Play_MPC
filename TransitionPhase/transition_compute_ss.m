
function [xs, us, alpha] = transition_compute_ss(x0, N, paramBefore, paramAfter, target)
    M = paramAfter.nb_subsystems;%length(param.activeDGU);
    %% create variables for optimizer
    nx = paramAfter.ni;
    nu = paramAfter.nu; % size 1
    % Input cell array of size N-1, each cell is an array of size nu*M
    U = sdpvar(repmat(nu,1,2*N-1), repmat(M,1,2*N-1),'full');
    % State cell array of size N (timestep), each cell is an array of size nx*M
    X = sdpvar(repmat(nx,1,2*N),repmat(M,1,2*N),'full'); % contains state of each subsystem i
    %X0 = sdpvar(nx,M,'full'); % state as rows and system number as column
    %X0 = sdpvar(repmat(nx,1,M),ones(1,M), 'full');
    X_Ni = cell(M,2*N-1); % cell array for neighbor states of i
    % Equilibrium state and input
    Us = sdpvar(nu, M,'full');
    Xs = sdpvar(nx, M,'full');
    X_sNi = cell(M,1);  % equilibrium neighbor state
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
    constraints = [];
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
    
    for i= union(paramBefore.activeDGU,paramAfter.activeDGU) % loop over active + to be plugged in/out
        constraints = [constraints, X{1}(:,i) == x0(:,i)];
        ni = size(paramBefore.A_Ni{i},1); 
        n_Ni = size(paramBefore.A_Ni{i},2); % get size of set of Neighbors
        % obtain sorted list of neighbors of system i
        neighbors_i = sort([i;neighbors(paramBefore.NetGraph, i)]);
        X_sNi{i} = sdpvar(n_Ni,1,'full'); % neighbor equilibrium state i
        %% Equilibrium constraints for system before plug in/plug out
        constraints = [constraints, Xs(:,i) == paramBefore.A_Ni{i}*X_sNi{i} + ...
                                                paramBefore.Bi{i}*Us(:,i)];
        % constraint for X_sNi =  concat of neighbor x_si
        constraints = [constraints, X_sNi{i} == ...
                                    reshape(Xs(:,neighbors_i),[],1)];
        constraints = [constraints, Us(:,i) == paramAfter.Ki{i}*Xs(:,i)+di(:,i)];
        constraints = [constraints, paramBefore.Gx_i{i}*Xs(:,i) <= paramBefore.fx_i{i}];
        %% Planning Horizon Loop
        for n = 1:N-1 
            % Neighbor States for each ith sytem at the kth horizon iteration
            X_Ni{i,n} = sdpvar(n_Ni,1,'full'); % neighbor set of state i
            % add a constraint for the neighbor state i to be equal to the
            % concatenated subsystem neighbor state vectors
            constraints = [constraints, X_Ni{i,n} == ...
                                        reshape(X{n}(:,neighbors_i),[],1)];
           
            % Distributed Dynamics
            constraints = [constraints, X{n+1}(:,i) == paramBefore.A_Ni{i}*X_Ni{i,n}+...
                                                       paramBefore.Bi{i}*U{n}(:,i)];
            % State and input constraints
            constraints = [constraints, paramBefore.Gx_i{i} * X{n}(:,i)...
                                      <= paramBefore.fx_i{i}];
            constraints = [constraints, paramBefore.Gu_i{i} * U{n}(:,i)...
                                       <= paramBefore.fu_i{i}];
            
        end
        % Objective
        if target == "reference"
            objective = objective + (Xs(:,i)-paramBefore.Xref{i})'*...
                                    (Xs(:,i)-paramBefore.Xref{i});
        elseif target == "current state"
            objective = objective + (Xs(:,i)-x0(:,i))'*...
                                    (Xs(:,i)-x0(:,i));
        else
            disp("objective not well defined, choose reference or current state");
        end
        % Terminal steady state condition
        constraints = [constraints, X{N}(:,i) == Xs(:,i)];    
    end
    
    for i=paramAfter.activeDGU % loop over active DGU after PnP operation
        %% Terminal Set conditions with reconfigurable terminal ingredients
        ni = size(paramAfter.A_Ni{i},1); 
        n_Ni = size(paramAfter.A_Ni{i},2); % get size of set of Neighbors
        % obtain sorted list of neighbors of system i
        neighbors_i = sort([i;neighbors(paramAfter.NetGraph, i)]);
        lambda{i} = sdpvar(n_Ni,1, 'full'); 
        %% Constraints corresponding to Eq 10, 11 and 12 in the Paper
        % Online Computation of Terminal Ingredients in Distributed Model
        % Predictive Control for Reference Tracking
        sumLambdaP_ij = 0;
        sumLambda_ij = 0;
        AbsSumLambdaP_ij = 0;
        % define alpha_Ni, c_Ni and Pij
        for j = 1:length(neighbors_i)
             alpha_Ni{i} = alpha_Ni{i} + ...
                            alpha(neighbors_i(j))*(paramAfter.Wij{i}{neighbors_i(j)})'*... 
                            (paramAfter.Wij{i}{neighbors_i(j)});
             c_Ni{i} = c_Ni{i} + paramAfter.Wij{i}{neighbors_i(j)}'*ci(:,neighbors_i(j));
             Pij = (paramAfter.Wij{i}{neighbors_i(j)})'*paramAfter.Pi{neighbors_i(j)}...
                    *(paramAfter.Wij{i}{neighbors_i(j)});
             sumLambdaP_ij = sumLambdaP_ij + ...
                             lambda{i}(j)*Pij;
             AbsSumLambdaP_ij = AbsSumLambdaP_ij + lambda{i}(j)*abs(Pij);
             sumLambda_ij = sumLambda_ij + lambda{i}(j);                   
        end
      
        % Equation 14: Approx of LMI with diagonal dominance
        PiInv = inv(paramAfter.Pi{i});
        constraints = [constraints, -bi(:,i) <= (paramAfter.A_Ni{i}+paramAfter.Bi{i}*paramAfter.K_Ni{i})...
                       *c_Ni{i} + paramAfter.Bi{i}*di(:,i) - ci(:,i) <= bi(:,i)];
        constraints = [constraints, alpha(i)-sumLambda_ij >= sum(bi(:,i))];    
        for k= 1:ni
            nondiag1 = sum(abs(PiInv(k,:))*alpha(i))- abs(PiInv(k,k))*alpha(i)+...
                       sum(abs(paramAfter.A_Ni{i}(k,:)+ paramAfter.Bi{i}(k)*paramAfter.K_Ni{i}(:))*alpha_Ni{i})...
                       + bi(k,i);
            constraints = [constraints, PiInv(k,k)*alpha(i) >= nondiag1];       
        end
        for k=1:n_Ni
            nondiag2 = sum(AbsSumLambdaP_ij(k,:)) - AbsSumLambdaP_ij(k,k) + ...
                       sum(abs(paramAfter.A_Ni{i}(:,k)+ paramAfter.Bi{i}(:)*paramAfter.K_Ni{i}(k))*alpha(i));
            constraints = [constraints, sumLambdaP_ij(k,k) >= nondiag2];
        end
    %   Equation 11
        for k=1:size(paramAfter.Gx_Ni{i},1) 
            sum_GxNorm2 = 0;
            for j=1:length(neighbors_i)
                sum_GxNorm2 = sum_GxNorm2 + norm(paramAfter.Gx_Ni{i}(k,:)*...
                             (paramAfter.Wij{i}{neighbors_i(j)})'*...
                              paramAfter.Pi{neighbors_i(j)}^(-1/2),2) * alpha(neighbors_i(j));
                
            end
            constraints = [constraints, paramAfter.Gx_Ni{i}(k,:)*c_Ni{i} + sum_GxNorm2 ...
                            <= paramAfter.fx_Ni{i}(k)];    
        end
       % Equation 12
        for k=1:size(paramAfter.Gu_i{i},1) 
            sum_GuNorm2 = 0;
            for j=1:length(neighbors_i)
                sum_GuNorm2 = sum_GuNorm2 + norm(paramAfter.Gu_i{i}(k,:)*paramAfter.K_Ni{i}*...
                             (paramAfter.Wij{i}{neighbors_i(j)})'*...
                              paramAfter.Pi{neighbors_i(j)}^(-1/2),2) * alpha(neighbors_i(j));
                
            end
            constraints = [constraints, paramAfter.Gu_i{i}(k,:)*paramAfter.K_Ni{i}*c_Ni{i} + ...
                           paramAfter.Gu_i{i}(k,:)*di(:,i) + sum_GuNorm2...
                           <= paramAfter.fu_i{i}(k)];    
        end
        %% Horizon Loop
        for n= N:2*N-1 % start from N here ? error in paper ?
            X_Ni{i,n} = sdpvar(n_Ni,1,'full'); % neighbor set of state i
            constraints = [constraints, X_Ni{i,n} == ...
                                        reshape(X{n}(:,neighbors_i),[],1)];
            % Distributed Dynamics
            constraints = [constraints, X{n+1}(:,i) == paramAfter.A_Ni{i}*X_Ni{i,n}+...
                                                       paramAfter.Bi{i}*U{n}(:,i)];
            constraints = [constraints, paramAfter.Gx_i{i} * X{n}(:,i)...
                                      <= paramAfter.fx_i{i}];
%             constraints = [constraints, paramAfter.Gx_Ni{i} * X_Ni{i,n}...
%                                   <= paramAfter.fx_Ni{i}];
           constraints = [constraints, paramAfter.Gu_i{i} * U{n}(:,i)...
                                       <= paramAfter.fu_i{i}];    
        end
        %%  Terminal Set condition
        constraints = [constraints, norm(paramAfter.Pi{i}^(1/2)*(X{2*N}(:,i)-ci(i)),2) <= alpha(i)];
          % constraints = [constraints, cone(paramAfter.Pi{i}^(1/2)*(X{2*N}(:,i)-ci(i)),alpha(i))];
        constraints = [constraints, alpha(i) >= 0];
    end  
    
    %% Create optimizer object 
    ops = sdpsettings('solver', 'MOSEK', 'verbose',2); %options
    parameters_in = [];
    solutions_out = {Xs, Us, alpha}; 
    optimizerObj = optimizer(constraints,objective,ops,parameters_in,solutions_out);
    %optim_results = optimizerObj(x0);
    [optim_results, ~, ~, ~, ~, feasibility] =  optimizerObj();
    if ~(feasibility.problem)
        disp("Feasible steady-state found");
    else
        error("Steady-state not found: P&P rejected");
    end
    xs = optim_results{1};
    us = optim_results{2};
    alpha = optim_results{3};
    
end