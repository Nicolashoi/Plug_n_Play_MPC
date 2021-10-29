function [Pi, P_Ni, K_Ni, Gamma_Ni] = terminal_costs_lyapunov_based(Q_Ni, Ri, system)
    objective = 0;
    % system choice
    if system == "COUPLED_OSCI"
        param = param_coupled_oscillator;
    else
        error("not implemented yet");
    end
    mi = param.size_subsystem;
    M = param.number_subsystem;
    %% decision variables
    global Ei E H_Ni Y_Ni E_Ni Ebar
    Ei = sdpvar(repmat(mi,1,M),repmat(mi,1,M)); % cell array of dimension M, 
    % each cell array is mi x mi 
    E = blkdiag(Ei{:}); 
    H_Ni = cell(1,M); Y_Ni = cell(1,M);
    
    %% dependent variables
    E_Ni = cell(1,M);
    Ebar = cell(1,M);
    
    constraints = define_constraints(Q_Ni, Ri, param, mi,M);
    %% OPTIMIZER
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0);
    %ops.mosek.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON';
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
        error('MOSEK solver thinks algorithm 1 is infeasible')
    end
    %% Map
    for i = M:-1:1
        Pi{i} = inv(value(Ei{i}));
        P_Ni{i} = inv(value(E_Ni{i}));
        K_Ni{i} = value(Y_Ni{i})*P_Ni{i};
        Gamma_Ni{i} = P_Ni{i}*value(H_Ni{i}) * P_Ni{i};
    end
end

function constraints = define_constraints(Q_Ni, Ri, param, mi, M)
    epsilon = 1e-5; % tolerance for positive definite constraint on E and E_Ni
    global Ei E H_Ni Y_Ni E_Ni Ebar
    A_Ni = cell(1,M);
    pi = size(param.Bi{1},2);
    constraints = [];
    %% Constraint definition
    for i = 1:M
        A_Ni{i} = param.U{i}*param.global_sysd.A*param.W{i};
        m_Ni = size(param.W{i},1); % size of the set Ni (i and it's neighbors)
        % decision variables dependent on it's subsystem neighbors set size
        Ebar{i} = param.W{i}*param.U{i}' * Ei{i}* param.U{i}* param.W{i}';
        E_Ni{i} = param.W{i}*E*param.W{i}';
        H_Ni{i} = sdpvar(m_Ni, m_Ni, 'full');
        Y_Ni{i} = sdpvar(pi, m_Ni, 'full');
        constraints = [constraints, Ei{i} >= epsilon*eye(size(Ei{i})),...
                       E_Ni{i} >= epsilon*eye(size(E_Ni{i}))] ; %P.D
        LMI{1} = [Ebar{i} + H_Ni{i}, E_Ni{i}*A_Ni{i}'+Y_Ni{i}'*param.Bi{i}',...
                 E_Ni{i}*Q_Ni{i}^(1/2), Y_Ni{i}'*Ri{i}^(1/2)];
        LMI{2} = [A_Ni{i}*E_Ni{i}+ param.Bi{i}*Y_Ni{i}, Ei{i}, zeros(mi,m_Ni),...
                  zeros(mi,pi)];
        LMI{3} = [(Q_Ni{i}^1/2)*E_Ni{i},zeros(m_Ni,mi), eye(m_Ni),...
                    zeros(m_Ni, pi)];
        LMI{4} = [(Ri{i}^1/2)*Y_Ni{i}, zeros(pi, mi), zeros(pi,m_Ni), eye(pi)];
        LMI_1 = [LMI{1}; LMI{2}; LMI{3}; LMI{4}];
        constraints = [constraints, LMI_1 >= 0];
        if i == 1
            LMI_2 = param.W{i}'*H_Ni{i}*param.W{i};
        else
            LMI_2 = LMI_2 + param.W{i}'*H_Ni{i}*param.W{i};
        end
    end
    constraints = [constraints, LMI_2 <= 0];
end