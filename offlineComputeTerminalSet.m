%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:
%  

function [Pi, Gamma_Ni, alpha_i] = offlineComputeTerminalSet(Q_Ni, Ri, param, passivity)
    % system choice
    M = length(param.activeDGU);%param.nb_subsystems;
    if ~passivity
        [Pi, P_Ni, K_Ni, Gamma_Ni] = terminal_costs_lyapunov_based(Q_Ni, Ri, param);
    elseif passivity
        Ki = cell(1,M); Pi = cell(1,M); Gamma_i = cell(1,M);
        K_Ni = cell(1,M); P_Ni = cell(1,M); Gamma_Ni = cell(1,M);
        for i=param.activeDGU
            [Ki{i}, ~, Pi{i}, Gamma_i{i}] = controller_passivity(param.Ai{i},...
                                            param.Bi{i}, param.Ci{i},param.Fi{i},...
                                            param.L_tilde, param.global_sysd.C, i);
            sprintf("passivity gain of system %d is", i)
            disp(Ki{i});
        end
        for i=param.activeDGU
            out_neighbors = sort([i; neighbors(param.NetGraph, i)]); 
            P_Ni{i} = blkdiag(Pi{out_neighbors}); 
            Gamma_Ni{i} = blkdiag(Gamma_i{out_neighbors});
            K_block = blkdiag(Ki{out_neighbors}); %block matrix of all neighbors K
            K_Ni{i} = K_block(out_neighbors==i,:); % extract only row corresponding to subsystem i
        end
    else
       error("ERROR: chose if passivity is to be used or not"); 
    end
    %% LP (equation 32 of the paper)
    constraints = []; % initialize constraints
    alpha = sdpvar(1,1,'full');
    for i = param.activeDGU

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
    alpha_i = value(alpha)/length(param.activeDGU);

end