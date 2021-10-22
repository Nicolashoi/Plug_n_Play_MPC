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


function [Pi, Gamma_i, alpha_i] = offline_distributed_MPC(Q_Ni, Ri, system)

if system == "COUPLED_OSCI"
    param = param_coupled_oscillator;
else
    error("not implemented yet");
end
B{1} = param.B_1;
B{2} = param.B_2;
constraints = [];
objective = 0;
mi = param.number_subsystem
pi = size(B{1},2);

for i = 1:mi
    A_Ni{i} = param.U{i]*param.global_sysd.A*param.W{i};
    m_Ni = size(param.W{i},1);
    Ebar{i} = sdpvar(m_Ni, m_Ni);
    E_Ni{i} = sdpvar(m_Ni, m_Ni);
    E{i} = sdpvar(mi,mi);
    H_Ni{i} = sdpvar(m_Ni, m_Ni, 'full');
    Y_Ni{i} = sdpvar(mi, m_Ni, 'full');
    LMI_1 = [Ebar{i} + H_Ni{i}, E_Ni{i}*A_Ni{i}'+Y_Ni{i}'*B{i}', E_Ni{i}*Q_Ni{i}^(1/2),...
           Y_Ni{i}'*Ri{i}^(1/2); A_Ni{i}*E_Ni{i}+ B{i}*Y_Ni{i}, E{i}, zeros(mi,m_Ni), zeros(mi,pi);...
           (Q_Ni{i}^1/2)*E_Ni{i},zeros(m_Ni,m_Ni), eye(m_Ni), zeros(m_Ni, pi);...
           (R_Ni{i}^1/2)*Y_Ni{i}, zeros(pi, mi), zeros(pi,m_Ni), eye(pi)];
    constraints = [constraints, LMI_1 >= 0]
    if i == 1
        LMI_2 = param.W{i}'*H_Ni{i}*param.W{i};
    else
        LMI_2 = LMI_2 + param.W{i}'*H_Ni{i}*param.W{i};
    end
end
constraints = [constraints, LMI_2 <= 0];


  %% OPTIMIZER
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0);
    %ops.mosek.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON';
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
        disp('MOSEK solver thinks it is infeasible')
    end
end