%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:
%  

function [P, Gamma_Ni, alpha_i] = offline_distributed_MPC(Q_Ni, Ri, system)

constraints = [];
objective = 0;

% system choice
if system == "COUPLED_OSCI"
    param = param_coupled_oscillator;
else
    error("not implemented yet");
end



mi = param.size_subsystem;
M = param.number_subsystem;
pi = size(param.Bi{1},2);
Ei = sdpvar(repmat(mi,1,M),repmat(mi,1,M)); % cell array of dimension M, each cell array is mi x mi 
E = blkdiag(Ei{:}); 
epsilon = 1e-5; % tolerance for positive definite constraint on E and E_Ni

for i = M:-1:1 % descending order so that Matlab allocates the necesssary free space
    A_Ni{i} = param.U{i}*param.global_sysd.A*param.W{i};
    m_Ni = size(param.W{i},1);
    %Ebar{i} = sdpvar(m_Ni, m_Ni);
    Ebar{i} = param.W{i}*param.U{i}' * Ei{i}* param.U{i}* param.W{i}';
    %E_N{i} = sdpvar(m_Ni, m_Ni);
    %Ei{i} = sdpvar(mi,mi);
    E_N{i} = param.W{i}*E*param.W{i}';
    H_N{i} = sdpvar(m_Ni, m_Ni, 'full');
    Y_N{i} = sdpvar(pi, m_Ni, 'full');
    constraints = [constraints,Ei{i} >= epsilon*eye(size(Ei{i})), E_N{i} >= epsilon*eye(size(E_N{i}))] ;
    LMI_1 = [Ebar{i} + H_N{i}, E_N{i}*A_Ni{i}'+Y_N{i}'*param.Bi{i}', E_N{i}*Q_Ni{i}^(1/2),...
           Y_N{i}'*Ri{i}^(1/2); A_Ni{i}*E_N{i}+ param.Bi{i}*Y_N{i}, Ei{i}, zeros(mi,m_Ni), zeros(mi,pi);...
           (Q_Ni{i}^1/2)*E_N{i},zeros(m_Ni,mi), eye(m_Ni), zeros(m_Ni, pi);...
           (Ri{i}^1/2)*Y_N{i}, zeros(pi, mi), zeros(pi,m_Ni), eye(pi)];
    constraints = [constraints, LMI_1 >= 0];
    if i == M
        LMI_2 = param.W{i}'*H_N{i}*param.W{i};
    else
        LMI_2 = LMI_2 + param.W{i}'*H_N{i}*param.W{i};
    end
end
constraints = [constraints, LMI_2 <= 0];

  %% OPTIMIZER
ops = sdpsettings('solver', 'MOSEK', 'verbose',0);
%ops.mosek.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON';
diagnostics = optimize(constraints, objective, ops);
if diagnostics.problem == 1
    error('MOSEK solver thinks algorithm 1 is infeasible')
end
%P = cell(1,M); % M subsystem, so M different Pi
for i = M:-1:1
    P{i} = inv(value(Ei{i}));
    P_N{i} = inv(value(E_N{i}));
    K_N{i} = value(Y_N{i})*P_N{i};
    Gamma_Ni{i} = P_N{i}*value(H_N{i}) * P_N{i};
end

% LP (equation 32 of the paper)
constraints = []; % reinitialize constraints
clear diagnostics objective ops
alpha = sdpvar(1,1,'full');
for i = M:-1:1
    
    for j= 1:size(param.Gx_i{i},1)
    constraints = [constraints,(norm(P{i}^(1/2)*param.Gx_i{i}(j,:)')^2)*alpha...
                   <= (param.fx_i{i}(j,:))^2];
    end
    for j=1:size(param.Gu_i{i},1)
    constraints = [constraints,(norm(P_N{i}^(1/2)*K_N{i}'*param.Gu_i{i}(j,:)')^2)...
                   *alpha <= (param.fu_i{i}(j,:))^2];
    end
end
objective = -alpha;
ops = sdpsettings('verbose', 0);
diagnostics = optimize(constraints, objective, ops);
if diagnostics.problem == 1
    error('Solver thinks algorithm 2 is infeasible')
end
alpha_i = value(alpha)/M;

end