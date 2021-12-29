%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:
%  
function [obj, Gamma_Ni, alpha_i] = offlineComputeTerminalSet(Q_Ni, Ri, obj)
    % system choice
    M = length(obj.activeDGU);%param.nb_subsystems;
    [Pi, P_Ni, K_Ni, Gamma_Ni] = terminal_costs_lyapunov_based(Q_Ni, Ri, obj);
    % Set controllers 
    obj.Pi = Pi;
    obj.K_Ni = K_Ni;
    %% LP (equation 32 of the paper)
    constraints = []; % initialize constraints
    alpha = sdpvar(1,1,'full');
    for i = obj.activeDGU
        for j= 1:size(obj.Gx_i{i},1)
        constraints = [constraints,(norm(Pi{i}^(1/2)*obj.Gx_i{i}(j,:)')^2)*alpha...
                       <= (obj.fx_i{i}(j,:))^2];
        end
        for j=1:size(obj.Gu_i{i},1)
        constraints = [constraints,(norm(P_Ni{i}^(1/2)*K_Ni{i}'*obj.Gu_i{i}(j,:)')^2)...
                       *alpha <= (obj.fu_i{i}(j,:))^2];
        end
    end
    objective = -alpha;
    ops = sdpsettings('verbose', 0);
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
        error('Solver thinks algorithm 2 is infeasible')
    end
    alpha_i = value(alpha)/length(obj.activeDGU);

end

function [Pi, P_Ni, K_Ni, Gamma_Ni] = terminal_costs_lyapunov_based(Q_Ni, Ri, param)
    ni = param.ni;
    M = param.nb_subsystems;%length(param.activeDGU);
    global Ei E H_Ni Y_Ni E_Ni Ebar
    %% decision variables
    Ei = sdpvar(repmat(ni,1,M),repmat(ni,1,M)); % cell array of dimension M, 
    % each cell array is ni x ni 
    H_Ni = cell(1,M); Y_Ni = cell(1,M);
    
    %% dependent variables
    E_Ni = cell(1,M);
    Ebar = cell(1,M);
    E = blkdiag(Ei{:}); 
    objective = 0;
    [constraints, objective] = define_constraints(objective, Q_Ni, Ri, param, ni);
    %% OPTIMIZER
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0);
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
        error('MOSEK solver thinks algorithm 1 is infeasible')
    end
    %% Map
    for i = param.activeDGU
        Pi{i} = inv(value(Ei{i}));
        P_Ni{i} = inv(value(E_Ni{i}));
        K_Ni{i} = value(Y_Ni{i})*P_Ni{i};
        sprintf("Local feedback control laws K_N%d", i)
        disp(K_Ni{i});
        Gamma_Ni{i} = P_Ni{i}*value(H_Ni{i}) * P_Ni{i};
    end
end

function [constraints, objective] = define_constraints(objective, Q_Ni, Ri, param, ni)
    epsilon = 1e-5; % tolerance for positive definite constraint on E and E_Ni
    global Ei E H_Ni Y_Ni E_Ni Ebar
    pi = size(param.Bi{1},2);
    constraints = [];
    %% Constraint definition
    LMI_2 = 0;
    for i = param.activeDGU
        n_Ni = size(param.W{i},1); % size of the set Ni (i and it's neighbors)
        % decision variables dependent on it's subsystem neighbors set size
        Ebar{i} = param.W{i}*param.U{i}' * Ei{i}* param.U{i}* param.W{i}';
        E_Ni{i} = param.W{i}*E*param.W{i}';
        H_Ni{i} = sdpvar(n_Ni, n_Ni, 'full');
        Y_Ni{i} = sdpvar(pi, n_Ni, 'full');
        
        constraints = [constraints, Ei{i} >= epsilon*eye(size(Ei{i})),...
                       E_Ni{i} >= epsilon*eye(size(E_Ni{i}))] ; %P.D
        LMI{1} = [Ebar{i} + H_Ni{i}, E_Ni{i}*param.A_Ni{i}'+Y_Ni{i}'*param.Bi{i}',...
                 E_Ni{i}*Q_Ni{i}^(1/2), Y_Ni{i}'*Ri{i}^(1/2)];
        LMI{2} = [param.A_Ni{i}*E_Ni{i}+ param.Bi{i}*Y_Ni{i}, Ei{i}, zeros(ni,n_Ni),...
                  zeros(ni,pi)];
        LMI{3} = [(Q_Ni{i}^1/2)*E_Ni{i},zeros(n_Ni,ni), eye(n_Ni),...
                    zeros(n_Ni, pi)];
        LMI{4} = [(Ri{i}^1/2)*Y_Ni{i}, zeros(pi, ni), zeros(pi,n_Ni), eye(pi)];
        % LMI 25 of the paper
        LMI_1 = [LMI{1}; LMI{2}; LMI{3}; LMI{4}];
        constraints = [constraints, LMI_1 >= 0];
        % LMI 26 of the paper
        LMI_2 = LMI_2 + param.W{i}'*H_Ni{i}*param.W{i};
        % objective
        objective = objective - trace(H_Ni{i});
    end
    constraints = [constraints, LMI_2 <= 0];
end