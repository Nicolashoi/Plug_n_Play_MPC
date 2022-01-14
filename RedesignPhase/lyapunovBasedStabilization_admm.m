%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF: NOT FULLLY IMPLEMENTED 
%%%%%%%%%%%%%% DEPRECATED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [K_Ni, P_Ni, Pi, Gamma_Ni] = lyapunovBasedStabilization_admm(Q_Ni,...
%                                         Ri, param)
%     rho = 0.25;
%     TMAX = 40;
%     Tk = 0; k = 2; l=1;
%     for i=param.activeDGU
%         neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
%         z_Ni{i,1}.E_Ni = zeros(length(neighbors_i)*param.ni);
%         y_Ni{i,1}.E_Ni = zeros(length(neighbors_i)*param.ni);
%         z_Ni{i,1}.H_Ni = cell(1,length(neighbors_i));
%         y_Ni{i,1}.H_Ni = cell(1,length(neighbors_i));
%         for j = 1:length(neighbors_i)
%             z_Ni{i,1}.H_Ni{j} = zeros(size(param.A_Ni{neighbors_i(j)},2));
%             y_Ni{i,1}.H_Ni{j} = zeros(size(param.A_Ni{neighbors_i(j)},2));
%         end
%         z_Ni{i,1}.Y_Ni = zeros(param.nu, length(neighbors_i)*param.ni);
%         y_Ni{i,1}.Y_Ni = zeros(param.nu, length(neighbors_i)*param.ni);
%             
%     end
%      while(Tk < TMAX)
%         for i=param.activeDGU % for each subsystems
%             % Solve a local optimizati% Loop while terminal terminal time not overrunnedon problem for each subsystem i
%             [w_Ni{i,k}, vi{i,k}, elapsedTime] = local_optim(i,k, Q_Ni{i}, Ri{i},...
%                                                 param, z_Ni{i,l}, y_Ni{i,l}, rho);
%             Tk = Tk + elapsedTime;
%         end
%      end
% end
% 
% function [w_Ni, vi, elapsedTime] = local_optim(i, k, Q_Ni, Ri, param, z_Ni, y_Ni, rho)
%     persistent localOptimizer
%     if k==2
%         localOptimizer{i} = init_optimizer(i, Q_Ni, Ri, param,rho);
%     end
%                                 
%     [solutionSet, ~, ~, ~, ~, optimTime] = localOptimizer{i}(z_Ni.E_Ni, z_Ni.H_Ni, z_Ni.Y_Ni, ...
%                      y_Ni.E_Ni, y_Ni.H_Ni, y_Ni.Y_Ni);
%     w_Ni.E_Ni = solutionSet{2};
%     w_Ni.H_Ni = solutionSet{3};
%     w_Ni.Y_Ni= solutionSet{4};
%     vi.Ei = solutionSet{1};    
%     elapsedTime= optimTime.solvertime;
% end
% 
% function localOptimizer = init_optimizer(i, Q_Ni, Ri, param, rho)
%     objective = 0;
%     constraints = [];
%     epsilon = 1e-5; % tolerance for positive definite constraint on E and E_Ni
%     neighbors_i = sort([i;neighbors(param.NetGraph, i)]);
%     % decision variables
%     n_Ni = size(param.A_Ni{i},2); % size of Neighbors set
%     Ei = sdpvar(param.ni); 
%     Ebar = param.W{i}*param.U{i}' * Ei* param.U{i}* param.W{i}';
%     E_Ni = sdpvar(n_Ni, n_Ni);
%     
% %     H_Ni = sdpvar(repmat(n_Ni,1, length(neighbors_i)), ...
% %                   repmat(n_Ni,1, length(neighbors_i)), 'full');
% 
%     Y_Ni = sdpvar(param.nu, n_Ni, 'full');
%     idx_i = neighbors_i == i;
%     % parameters that will change as inputs
%     z_Ni.E_Ni = sdpvar(n_Ni, n_Ni);
%     y_Ni.E_Ni = sdpvar(n_Ni, n_Ni);
%     H_Ni = cell(1,length(neighbors_i));
%     z_Ni.H_Ni = cell(1,length(neighbors_i));
%     y_Ni.H_Ni = cell(1,length(neighbors_i));
%     for j = 1:length(neighbors_i)
%         n_Nj = size(param.A_Ni{neighbors_i(j)},2);
%         H_Ni{j} = sdpvar(n_Nj, n_Nj, 'full');
%         z_Ni.H_Ni{j} = sdpvar(n_Nj, n_Nj, 'full');
%         y_Ni.H_Ni{j} = sdpvar(n_Nj, n_Nj, 'full');
%     end
%    
%     z_Ni.Y_Ni = sdpvar(param.nu, n_Ni, 'full');
%     y_Ni.Y_Ni = sdpvar(param.nu, n_Ni, 'full');
%     
%     
%     % constraints
%     constraints = [constraints, Ei >= epsilon*eye(size(Ei)),...
%                        E_Ni >= epsilon*eye(size(E_Ni))] ; %P.D
%         
%     LMI{1} = [Ebar + H_Ni{idx_i}, E_Ni*param.A_Ni{i}'+Y_Ni'*param.Bi{i}',...
%                  E_Ni*Q_Ni^(1/2), Y_Ni'*Ri^(1/2)];
%     LMI{2} = [param.A_Ni{i}*E_Ni+ param.Bi{i}*Y_Ni, Ei, zeros(param.ni,n_Ni),...
%               zeros(param.ni,param.nu)];
%     LMI{3} = [Q_Ni^1/2*E_Ni,zeros(n_Ni,param.ni), eye(n_Ni),...
%                 zeros(n_Ni, param.nu)];
%     LMI{4} = [Ri^1/2*Y_Ni, zeros(param.nu, param.ni), zeros(param.nu,n_Ni),...
%                eye(param.nu)];
%     LMI_1 = vertcat(LMI{:});
%         %LMI_1 = [LMI{1}; LMI{2}; LMI{3}; LMI{4}];
%     constraints = [constraints, LMI_1 >= 0];
%    
%     LMI_2 = 0;
%     for j = neighbors_i'
%         idx_j = neighbors_i==j;
%         LMI_2 = LMI_2 + param.W{j}'*H_Ni{j}*param.W{j};
%         objective = objective + sum(y_Ni.H_Ni{j}'*(H_Ni{j} - z_Ni.H_Ni{j}) + ...
%                     rho/2 * (H_Ni{j} - z_Ni.H_Ni{j})'*(H_Ni{j} - z_Ni.H_Ni{j}), 'all');
%     end
%     constraints = [constraints, LMI_2 <= 0];
%     
%     objective = objective + sum(y_Ni.E_Ni'*(E_Ni-z_Ni.E_Ni) + y_Ni.Y_Ni'*...
%                 (Y_Ni-z_Ni.Y_Ni) + rho/2 *(E_Ni-z_Ni.E_Ni)'*(E_Ni-z_Ni.E_Ni)+...
%                  rho/2 *(Y_Ni-z_Ni.Y_Ni)'*(Y_Ni-z_Ni.Y_Ni), 'all');
%     parameters_in = {z_Ni.E_Ni, z_Ni.H_Ni{1}, z_Ni.H_Ni{2}, z_Ni.Y_Ni, ...
%                      y_Ni.E_Ni, y_Ni.H_Ni{1},y_Ni.H_Ni{2}, y_Ni.Y_Ni};
%     ops = sdpsettings('solver', 'MOSEK', 'verbose',1); %options
%     solutions_out = {Ei, E_Ni, [H_Ni{:}], Y_Ni};
%     optimizer(constraints,objective,ops,parameters_in,solutions_out);
% end
% 
% 
% 
% 
% function zi = update_global_copy(wi)
%     fn = fieldnames(wi{1});
%     for ii = 1:numel(fn)
%         % extract fieldname structure for every neighbors
%         extract = cellfun(@(x) x.(fn{ii}), wi(:), 'Un', false); 
%         % concatenate the same structure fieldname into 3rd dimension
%         zi.(fn{ii}) = mean(cat(3, extract{:}),3); % mean along 3rd dimension
%     end
% end
% 
% function result = diff_struct(w_Ni, z_Ni)
%     fn = fieldnames(w_Ni);
%     for i = 1:numel(fn)
%         result.(fn{i}) = w_Ni.(fn{i}) - z_Ni.(fn{i});
%     end
% end
% 
% function result = add_struct(struct1, struct2)
%     fn = fieldnames(struct1);
%     for i = 1:numel(fn)
%         result.(fn{i}) = struct1.(fn{i}) + struct2.(fn{i});
%     end
% end