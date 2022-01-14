function [Qi, Ri, decVariables] = computeQi_Ri(param, i)
    constraints = [];
    Ri = sdpvar(param.nu, param.nu);
    Qi = sdpvar(param.ni, param.ni); % symmetric positive definite
    zDec = sdpvar(param.nb_subsystems,1, 'full');
    neighbors_i = sort([i; neighbors(param.NetGraph, i)]);
    Plift_cell = cell(length(neighbors_i),1); 
    Acl = param.A_Ni{i} + param.Bi{i}*param.K_Ni{i};
    constraints = [constraints, zDec(neighbors_i) >= 0];
 
    for j = neighbors_i'
        Plift_cell{j} = param.Pi{j}*zDec(j);
    end
    PiLift = blkdiag(Plift_cell{:});
    idx_i_mat = diag(neighbors_i == i);
    QiLift = kron(idx_i_mat, Qi);
    
    LMI_i = PiLift - QiLift - Acl'*param.Pi{i}*Acl - ...
            param.K_Ni{i}'*Ri*param.K_Ni{i};
    
    constraints = [constraints, zDec(setdiff(1:param.nb_subsystems, neighbors_i))...
                                == 0];

    constraints = [constraints, Qi >= 1e-2*eye(param.ni), ...
                                Ri >= 1e-2*eye(param.nu)];
    constraints = [constraints, LMI_i>=0];
    
    objective = sum(zDec(neighbors_i));
    objective = objective - 0.1*trace(Qi) - 0.1*trace(Ri);
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0); %options
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
        sprintf('MOSEK solver thinks it is infeasible for system %d', i)
    end
    
    Qi = value(Qi);
    Ri = value(Ri);
    decVariables = value(zDec);
end
