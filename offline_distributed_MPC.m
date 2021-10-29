%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:
%  

function [Pi, Gamma_Ni, alpha_i] = offline_distributed_MPC(Q_Ni, Ri, system)
    % system choice
    if system == "COUPLED_OSCI"
            param = param_coupled_oscillator;
        else
            error("not implemented yet");
    end
    
    M = param.number_subsystem;

    [Pi, P_Ni, K_Ni, Gamma_Ni] = terminal_costs_lyapunov_based(Q_Ni, Ri, system);

    %% LP (equation 32 of the paper)
    constraints = []; % initialize constraints
    alpha = sdpvar(1,1,'full');
    for i = M:-1:1

        for j= 1:size(param.Gx_i{i},1)
        constraints = [constraints,(norm(Pi{i}^(1/2)*param.Gx_i{i}(j,:)')^2)*alpha...
                       <= (param.fx_i{i}(j,:))^2];
        end
        for j=1:size(param.Gu_i{i},1)
        constraints = [constraints,(norm(P_Ni{i}^(1/2)*K_Ni{i}'*param.Gu_i{i}(j,:)')^2)...
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