function [Q_Ni, Ri] = computeQi_Ri(param, i)
    R = sdpvar(param.nu, param.nu);
    Qi = sdpvar(param.ni, param.ni); % symmetric positive definite
    epsilon = 1e-5;
    %Qi_lift=sdpvar(size(param.A_Ni{i},2));
    out_neighbors = sort([i; neighbors(param.NetGraph, i)]); 
    idx_i_mat = diag(out_neighbors == i);
    Qi_lift = kron(idx_i_mat, Qi);
    Pi_lift = kron(idx_i_mat, param.Pi{i});
    LHS = (param.A_Ni{i}+ param.Bi{i}*param.K_Ni{i})'*param.Pi{i} *...
          (param.A_Ni{i}+ param.Bi{i}*param.K_Ni{i}) - Pi_lift;
    RHS = -Qi_lift - param.K_Ni{i}'*R*param.K_Ni{i};
    
    constraints = [Qi >= epsilon*eye(param.ni),R >= epsilon*eye(param.nu),...
                   LHS <= RHS];
    objective = trace(Qi) + trace(R);
    
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0); %options
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
       sprintf('MOSEK solver thinks it is infeasible for system %d', i)
    end
    
    Q_Ni = value(Qi_lift);
    Ri = value(R);
end