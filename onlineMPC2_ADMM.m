function u0 = onlineMPC2_ADMM(x0,Q_Ni, Ri, N, param)
    p = 1/2;
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
        tstart  =tic;
    %% Constraints: Outer loop over subsystems, inner loop over Horizon
        for i=param.activeDGU % loop over all subsystems
            objective_i = 0;
            Xi = sdpvar(param.ni,N, 'full');
            Ui = sdpvar(param.nu,N-1, 'full');
            n_Ni = size(param.A_Ni{i},2); % get size of set of Neighbors
            X_eNi = sdpvar(n_Ni,1,'full'); % neighbor equilibrium state i
            X_Ni = sdpvar(n_Ni, N, 'full');
            Xei = sdpvar(ni,1,'full');
            Uei = sdpvar(nu,1,'full');
            di = sdpvar(nu,1,'full');
            constraints_i = [constraints_i, Xi(:,1) == x0(:,i)];
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
                objective_i = objective_i + ...
                            (X_Ni(:,n)-X_eNi)'*Q_Ni{i}*(X_Ni(:,n)-X_eNi)+...
                            (Ui(:,n)-Uei)'*Ri{i}*(Ui(:,n)-Uei) + ...
                                            
            end
            %% Terminal cost
%             objective = objective + ... 
%                         (Xi{end}(:,i)-Xe(:,i))'*Pi{i}*(Xi{end}(:,i)-Xe(:,i))+...
%                         (Xe(:,i) - param.Xref{i})'*Si*(Xe(:,i) - param.Xref{i});

            %%  Terminal Set condition
%             constraints_i = [constraints_i, (Xi{end}(:,i)-ci(i))'*Pi{i}*(Xi{end}(:,i)-ci(i))...
%                                         <= alpha(i)^2];      
          
            w_Ni{i,k}.x_Ni = X_Ni;
            w_Ni{i,k}.x_eNi = X_eNi;
            %wNi{i,k}.alpha_Ni = value(X_eNi);
            %wNi{i,k}.cNi = value(cNi);
            vi{i,k}.ui = Ui;
            vi{i,k}.uei = Uei;
            vi{i,k}.di = di;
            %vi{i,k}.lambda_ij = value(lambda_ij);
            wi{i,k,i}.xi = Xi;
            wi{i,k,i}.xei = Xei;
            %wi{i,k,i}.alpha_i = value(alpha_i); % ith value computed by system i
            %wi{i,k,i}.ci = value(ci);
            [w_Ni{i,k}, vi{i,k}] = lagrangian(N, p, constraints_i, objective_i, ...
                                   w_Ni{i,k}, vi{i,k}, z_Ni{i,l}, y_Ni{i,l});
                               
            for j = neighbors_i
                wi{i,k,j}.xi = param.Wij{i}{j}*w_Ni{i,k}.x_Ni(:,:);
                wi{i,k,j}.xei = param.Wij{i}{j}*w_Ni{i,k}.x_eNi;
            end
        end  
        for i = param.activeDGU
            zi{i,k} = update_global_copy(wi, i,k);
        end
        for i=param.activeDGU
            neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
            z_Ni{i,k}.x_Ni = vertcat(zi{neighbors_i,k}.xi ,[],1);
            z_Ni{i,k}.x_eNi = vertcat(zi{neighbors_i,k}.xei ,[],1);
            %z_Ni{i,k}.alpha_Ni = vertcat(z{neighbors_i,k}.alpha_i ,[],1);
            %z_Ni{i,k}.c_Ni = vertcat(z{neighbors_i,k}.ci ,[],1);
            
            y_Ni{i,k} = add_struct(y_Ni{i,l}, p.*diff_struct(w_Ni{i,k}, z_Ni{i,k}));  
        end
        Tk = Tk + toc(tStart);
    end
    
end

function zi = update_global_copy(wi, i, k}
    fn = fieldname(wi{i,k,i});
    for ii = 1:numel(fn)
        zi.(fn{ii}) = mean(horzcat(wi{i,k,:}.(fn(ii))),2);
%             z{i,k}.xi = mean(horzcat(wi{i,k,:}.xi),2);
%             z{i,k}.xei = mean(horzcat(wi{i,k,:}.xei),2); 
    end
end

function [w_Ni_new, vi_new] = lagrangian(p, constraints_i, objective_i, w_Ni, vi, z_Ni, y_Ni)
    %fn = fieldname(w_Ni);
    wNiMinuszNi = diff_struct(w_Ni, z_Ni);
    fn = fieldname(wNiMinuszNi);
    objective = objective_i;
    for i = 1:numel(fn)
        objective = objective + y_Ni.(fn{i}).*wNiMinuszNi.(fn{i}) + ...
                                p/2.* (wNiMinuszNi.(fn{i})).^2;
%         objective = objective + y_Ni.(fn{i}).*(w_Ni.(fn{i}) - z_Ni.(fn{i})) + ...
%                     p/2*(w_Ni.(fn{i}) - z_Ni.(fn{i})).^2;
    end
    ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
    diagnostics = optimize(constraints_i, objective, ops);
    if diagnostics.problem == 1
       fprintf("MOSEK solver thinks it is infeasible");
    end
    for i = 1:numel(fn)
        w_Ni_new.(fn{i}) = value(w_Ni.(fn{i}));
    end
    fnv = fieldname(vi);
    for ii = 1:numel(fnv)
        vi_new.(fnv{ii}) = value(vi.(fnv{ii}));
    end  
end

function result = diff_struct(w_Ni, z_Ni)
    fn = fieldname(w_Ni);
    for i = 1:numel(fn)
        result.(fn{i}) = w_Ni.(fn{i}) - z_Ni.(fn{i});
    end
end

function result = add_struct(struct1, struct2)
    fn = fieldname(struct1);
    for i = 1:numel(fn)
        result.(fn{i}) = struct1.(fn{i}) + struct2.(fn{i});
    end
end