%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:
%  
% INPUT:
%   Ai, Bi, Ci, Fi: LTI dynamics of subsystem i
%   U: L_tilde*C
%   W: C^T * L_tilde where L_tilde is the augmented laplacian
% OUTPUT:
%   Ki: decentralized passivity-based control state feedback gain
%   Di: positive-definite diagonal matrix
%   Pi: positive-definite matrix
%   Gamma_i: positive-definite diagonal matrix


function [P, Gamma_Ni, alpha_i] = offline_distributed_MPC(Q_Ni, Ri, system)

if system == "COUPLED_OSCI"
    param = param_coupled_oscillator;
else
    error("not implemented yet");
end
B{1} = param.B_1;
B{2} = param.B_2;
constraints = [];
objective = 0;
mi = param.size_subsystem;
M = param.number_subsystem;
pi = size(B{1},2);
E = sdpvar(repmat(mi,1,mi),repmat(mi,1,mi));
%Ebar = sdpvar(repmat(m_Ni,1,mi),repmat(m_Ni,1,mi)); cell array of mi
%matrices of dimension m_Ni x m_Ni
epsilon = 1e-5;
for i = 1:M
    A_Ni{i} = param.U{i}*param.global_sysd.A*param.W{i};
    m_Ni = size(param.W{i},1);
    Ebar{i} = sdpvar(m_Ni, m_Ni);
    E_N{i} = sdpvar(m_Ni, m_Ni);
    E{i} = sdpvar(mi,mi);
    H_N{i} = sdpvar(m_Ni, m_Ni, 'full');
    Y_N{i} = sdpvar(pi, m_Ni, 'full');
    constraints = [constraints,E{i} >= epsilon*eye(size(E{i})), E_N{i} >= epsilon*eye(size(E_N{i}))] 
    LMI_1 = [Ebar{i} + H_N{i}, E_N{i}*A_Ni{i}'+Y_N{i}'*B{i}', E_N{i}*Q_Ni{i}^(1/2),...
           Y_N{i}'*Ri{i}^(1/2); A_Ni{i}*E_N{i}+ B{i}*Y_N{i}, E{i}, zeros(mi,m_Ni), zeros(mi,pi);...
           (Q_Ni{i}^1/2)*E_N{i},zeros(m_Ni,mi), eye(m_Ni), zeros(m_Ni, pi);...
           (Ri{i}^1/2)*Y_N{i}, zeros(pi, mi), zeros(pi,m_Ni), eye(pi)];
    constraints = [constraints, LMI_1 >= 0];
    if i == 1
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
P = cell(mi,mi,M); % M subsystem, so M different Pi
for i = 1:M
    P{i} = inv(value(E{i}));
    P_N{i} = inv(value(E_N{i}));
    K_N{i} = value(Y_N{i})*P_N{i};
    Gamma_Ni{i} = P_N{i}*value(H_N{i}) * P_N{i};
end

% LP (equation 32 of the paper)
constraints = []; % reinitialize constraints
clear diagnostics objective ops
alpha = sdpvar(1,1,'full');
for i = 1:M
    for j= 1:size(param.Gx{i},1)
    constraints = [constraints,(norm(P{i}^(1/2)*param.Gx{i}(j,:)')^2)*alpha...
                   <= (param.fx{i}(j,:))^2];
    end
    for j=1:size(param.Gu{i},1)
    constraints = [constraints,(norm(P_N{i}^(1/2)*K_N{i}'*param.Gu{i}(j,:)')^2)...
                   *alpha <= (param.fu{i}(j,:))^2];
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