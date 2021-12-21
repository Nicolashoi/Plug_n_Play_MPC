function [Q, R] = computeQi_Ri_global(param)
    constraints = [];
    %R = sdpvar(param.nu*param.nb_subsystems, param.nu*param.nb_subsystems);
    R = diag(sdpvar(param.nb_subsystems,1));
    %Q = sdpvar(param.ni*param.nb_subsystems); % symmetric positive definite
    %Q = diag(sdpvar(param
    %A = param.global_sysd.A; 
    % B = param.global_sysd.B;
    B = blkdiag(param.Bi{:});
    A = [];
    for j = 1:param.nb_subsystems
        Qi{j} = sdpvar(param.ni);
        A = [A; param.A_Ni{j}*param.W{j}]; % concatenate horizontally
    end
    Q = blkdiag(Qi{:});
   
    Pblk = blkdiag(param.Pi{:});
    Kblk = blkdiag(param.Ki{:});
    LMI = Pblk - (A+B*Kblk)'*Pblk*(A+B*Kblk) - Q - Kblk'*R*Kblk;
    %constraints = [constraints, LMI >= 0, Q>= 1e-5*eye(12), R>= 1e-5*eye(6)];
    constraints = [constraints, Q>= 1e-5*eye(12), R>= 1e-5*eye(6)];
    for k=1:size(LMI,1)
        constraints = [constraints, LMI(k,k) >= ...
                      sum(abs(LMI(k,:)))-abs(LMI(k, k))];    
    end
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0); %options
    diagnostics = optimize(constraints, [], ops);
    if diagnostics.problem == 1
        sprintf('MOSEK solver thinks it is infeasible')
    end
    
    Q = value(Q);
    R = value(R);

end