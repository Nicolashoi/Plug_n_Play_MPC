function [Qi, Ri] = computeQi_Ri_old(param, i)
    constraints = [];
    Ri = sdpvar(param.nu, param.nu);
    Qi = sdpvar(param.ni, param.ni); % symmetric positive definite
    neighbors_i = sort([i; neighbors(param.NetGraph, i)]); 

    % At = [];
   % Aii = param.A_Ni{i}*param.Wij{i}{i}';
%     diagTerm = Aii'*param.Pi{i}*param.Bi{i}*param.Ki{i} + (param.Bi{i}*param.Ki{i})'...
%                 *param.Pi{i}*Aii + (param.Bi{i}*param.Ki{i})'*param.Pi{i}*...
%                 (param.Bi{i}*param.Ki{i}) + param.Pi{i};
    A = [];
    BK = [];
    for j = neighbors_i'
        Aij = param.A_Ni{j}*param.Wij{j}{i}'; % obtain Aij for system i from A_Ni
        %diagTerm = diagTerm + Aij'*param.Pi{j}*Aij; % add diagonal element of A^T*P*A
        A = [A; param.A_Ni{j}*param.W{j}]; % concatenate horizontally
        BK = [BK;param.Bi{j}*param.Ki{j}*param.U{j}]; 
%         if j==i
%             % if it is system i
%             BK = [BK;param.Bi{i}*param.Ki{i}*param.U{i}]; 
%         else
%             % For the neighbors where the row is not equal to system i we add zeros
%             BK = [BK; zeros(param.ni, param.nb_subsystems*param.ni)]; % padd zeros
%         end
%         At = [At Aij']; 
%         AtPBK = [AtPBK Aij'*param.Pi{i}*param.Bi{i}*param.Ki{i}];
    end
    Pblk = blkdiag(param.Pi{neighbors_i});
    Acl = (A+BK)'*Pblk*(A+BK);
    Acl = param.U{i}*Acl; % extract only rows for system i
    PiLift = param.Pi{i}*param.U{i};
    QiLift = Qi*param.U{i};
    Rlift = (param.Ki{i}'*Ri*param.Ki{i})*param.U{i};
    LMI_i = PiLift -QiLift - Rlift - Acl;
    shiftDiag = 0;
    for k=1:size(LMI_i,1)
        constraints = [constraints, LMI_i(k, i*param.ni-1+shiftDiag) >= ...
                      sum(abs(LMI_i(k,:)))-abs(LMI_i(k, i*param.ni-1+shiftDiag))];
        shiftDiag = shiftDiag+1;
    end
    constraints = [constraints, Qi >= 1e-2*eye(param.ni)];
    constraints = [constraints, Ri >= 1e-2 * eye(param.nu)];
   
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0); %options
    diagnostics = optimize(constraints, [], ops);
    if diagnostics.problem == 1
        sprintf('MOSEK solver thinks it is infeasible for system %d', i)
    end
    
    Qi = value(Qi);
    Ri = value(Ri);

end