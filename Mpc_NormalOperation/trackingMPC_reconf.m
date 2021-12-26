function [u0, alpha] = trackingMPC_reconf(x0, Q_Ni, Ri, N, param)
    persistent mpc_optimizer
    % initialize controller, if not done already
    if isempty(mpc_optimizer)
        mpc_optimizer = init_optimizer(Q_Ni, Ri, N, param);
    end
    [sol, ~, ~, ~, ~, feasibility]= mpc_optimizer(x0);
    if isequal(feasibility.problem,1)
        disp(feasibility.infostr);
    end
    u0 = sol{1};
    alpha = sol{2};
    c = sol{3};
    Xend = sol{4};

end
 
function mpc_optimizer = init_optimizer(Q_Ni, Ri, N, param)
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
    constraints = [];
    S = cell(M,1);
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
    
    for i=param.activeDGU % loop over all subsystems
        constraints = [constraints, X{1}(:,i) == X0(:,i)];
        S{i} = 1000*eye(size(param.Ai{i},1));
        ni = size(param.A_Ni{i},1); 
        n_Ni = size(param.A_Ni{i},2); % get size of set of Neighbors
        % obtain sorted list of neighbors of system i
        neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
        X_eNi{i} = sdpvar(n_Ni,1,'full'); % neighbor equilibrium state i
        lambda{i} = sdpvar(n_Ni,1, 'full');
        % constraint for X_eNi =  concat of neighbor x_ei
        constraints = [constraints, X_eNi{i} == ...
                                    reshape(Xe(:,neighbors_i),[],1)];
        %% Equilibrium constraints
        constraints = [constraints, Xe(:,i) == param.A_Ni{i}*X_eNi{i} + ...
                                                param.Bi{i}*Ue(:,i)];
        constraints = [constraints, Ue(:,i) == param.K_Ni{i}*X_eNi{i} + di(:,i)];  
        
        %% Constraints corresponding to Eq 10, 11 and 12 in the Paper
        % Online Computation of Terminal Ingredients in Distributed Model
        % Predictive Control for Reference Tracking
        sumLambdaP_ij = 0;
        sumLambda_ij = 0;
        AbsSumLambdaP_ij = 0;
        % define alpha_Ni, c_Ni and Pij
        for j = 1:length(neighbors_i)
%             constraints = [constraints,  alpha_Ni{i} == ...
%             alpha_Ni{i} + alpha(neighbors(j))*param.Wij{i}{neighbors(j)}*... 
%                             (param.Wij{i}{neighbors(j)})'];
%             constraints = [constraints, c_Ni{i} == ...  
%             c_Ni{i} + param.Wij{i}{neighbors(j)}'*ci(:,neighbors(j))];
             alpha_Ni{i} = alpha_Ni{i} + ...
                            alpha(neighbors_i(j))*(param.Wij{i}{neighbors_i(j)})'*... 
                            (param.Wij{i}{neighbors_i(j)});
             c_Ni{i} = c_Ni{i} + param.Wij{i}{neighbors_i(j)}'*ci(:,neighbors_i(j));
             Pij = (param.Wij{i}{neighbors_i(j)})'*param.Pi{neighbors_i(j)}...
                    *(param.Wij{i}{neighbors_i(j)});
             sumLambdaP_ij = sumLambdaP_ij + ...
                             lambda{i}(j)*Pij;
             AbsSumLambdaP_ij = AbsSumLambdaP_ij + lambda{i}(j)*abs(Pij);
             sumLambda_ij = sumLambda_ij + lambda{i}(j);                   
        end
        %% Equation 10
%         LMI_1 = [inv(param.Pi{i})*alpha(i), (param.A_Ni{i}+param.Bi{i}*param.K_Ni{i})*alpha_Ni{i},...
%                (param.A_Ni{i}+param.Bi{i}*param.K_Ni{i})*c_Ni{i} + param.Bi{i}*di(:,i) - ci(:,i)];
%         LMI_2 = [((param.A_Ni{i}+param.Bi{i}*param.K_Ni{i})*alpha_Ni{i})', sumLambdaP_ij,...
%                 zeros(size(sumLambdaP_ij,1), size(ci(:,i),2))];
%         LMI_3 = [((param.A_Ni{i}+param.Bi{i}*param.K_Ni{i})*c_Ni{i} + param.Bi{i}*di(:,i) - ci(:,i))',...
%                     zeros(size(alpha(i),1),size(sumLambdaP_ij,2)) , alpha(i) - sumLambda_ij];
%                  
%         LMI = [LMI_1; LMI_2; LMI_3];
%         constraints = [constraints, LMI>=0];    
        
        %% Equation 14: Approx of LMI with diagonal dominance
%         PiInv = inv(param.Pi{i});
%         constraints = [constraints, (param.A_Ni{i}+param.Bi{i}*param.K_Ni{i})...
%                        *c_Ni{i} + param.Bi{i}*di(:,i) - ci(:,i) == bi(:,i)];
%         constraints = [constraints, alpha(i)-sumLambda_ij >= sum(bi(:,i))];    
%         for k= 1:ni
%             nondiag1 = sum(abs(PiInv(k,:))*alpha(i))- abs(PiInv(k,k))*alpha(i)+...
%                        sum(abs(param.A_Ni{i}(k,:)+ param.Bi{i}(k)*param.K_Ni{i}(:))*alpha_Ni{i})...
%                        + bi(k,i);
%             constraints = [constraints, PiInv(k,k)*alpha(i) >= nondiag1];       
%         end
%         for k=1:n_Ni
%             nondiag2 = sum(AbsSumLambdaP_ij(k,:)) - AbsSumLambdaP_ij(k,k) + ...
%                        sum(abs(param.A_Ni{i}(:,k)+ param.Bi{i}(:)*param.K_Ni{i}(k))*alpha(i));
%             constraints = [constraints, sumLambdaP_ij(k,k) >= nondiag2];
%         end
        %% Equation 11
%         for k=1:size(param.Gx_Ni{i},1) 
%             sum_GxNorm2 = 0;
%             for j=1:length(neighbors_i)
%                 sum_GxNorm2 = sum_GxNorm2 + norm(param.Gx_Ni{i}(k,:)*...
%                              (param.Wij{i}{neighbors_i(j)})'*...
%                               param.Pi{neighbors_i(j)}^(-1/2),2) * alpha(neighbors_i(j));
%                 
%             end
%             constraints = [constraints, param.Gx_Ni{i}(k,:)*c_Ni{i} + sum_GxNorm2 ...
%                             <= param.fx_Ni{i}(k)];    
%         end
        %% Equation 12
%         for k=1:size(param.Gu_i{i},1) 
%             sum_GuNorm2 = 0;
%             for j=1:length(neighbors_i)
%                 sum_GuNorm2 = sum_GuNorm2 + norm(param.Gu_i{i}(k,:)*param.K_Ni{i}*...
%                              (param.Wij{i}{neighbors_i(j)})'*...
%                               param.Pi{neighbors_i(j)}^(-1/2),2) * alpha(neighbors_i(j));
%                 
%             end
%             constraints = [constraints, param.Gu_i{i}(k,:)*param.K_Ni{i}*c_Ni{i} + ...
%                            param.Gu_i{i}(k,:)*di(:,i) + sum_GuNorm2...
%                            <= param.fu_i{i}(k)];    
%         end
        
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
            constraints = [constraints, param.Gx_Ni{i} * X_Ni{i,n}...
                                      <= param.fx_Ni{i}];
%             constraints = [constraints, param.Gx_i{i} * X{n}(:,i)...
%                                       <= param.fx_i{i}];
            constraints = [constraints, param.Gu_i{i} * U{n}(:,i)...
                                       <= param.fu_i{i}];
            % Objective
            objective = objective + ...
                        (X_Ni{i,n}-X_eNi{i})'*Q_Ni{i}*(X_Ni{i,n}-X_eNi{i})+...
                        (U{n}(:,i)-Ue(:,i))'*Ri{i}*(U{n}(:,i)-Ue(:,i));                    
        end
        %% Terminal cost
        objective = objective + ... 
                    (X{end}(:,i)-Xe(:,i))'*param.Pi{i}*(X{end}(:,i)-Xe(:,i))+...
                    (Xe(:,i) - param.Xref{i})'*S{i}*(Xe(:,i) - param.Xref{i});
                        
        %%  Terminal Set condition
%         constraints = [constraints, (X{N}(:,i)-ci(:,i))'*param.Pi{i}*(X{N}(:,i)-ci(:,i))...
%                                     <= 2];   
        LMI_terminal = [inv(param.Pi{i})*alpha(i), X{end}(:,i) - ci(:,i);...
                        X{end}(:,i)' - ci(:,i)', alpha(i)];
         constraints = [constraints, LMI_terminal >= 0];
         
    end        
    %% Create optimizer object 
    ops = sdpsettings('solver', 'MOSEK', 'verbose',2); %options
    parameters_in = {X0};
    %solutions_out = {[U{:}], [X_eNi{1}], [X_eNi{2}], [X_eNi{3}], di, Ue}; % general form 
    solutions_out = {U{1}, alpha, ci, X{end}}; % get U0 for each subsystem, size nu x M
    mpc_optimizer = optimizer(constraints,objective,ops,parameters_in,solutions_out);
end