function [Qi, Ri, decVariables] = computeQi_Ri_new(param, i)
    constraints = [];
    Ri = sdpvar(param.nu, param.nu);
    Qi = sdpvar(param.ni, param.ni); % symmetric positive definite
    neighbors_i = sort([i; neighbors(param.NetGraph, i)]); 
    Pblk = blkdiag(param.Pi{neighbors_i});
    %Kblk = blkdiag(param.Ki{neighbors_i});
    BiLift = param.Wij{i}{i}'*param.Bi{i};
    zDec = sdpvar(param.nb_subsystems,1, 'full');
    Plift_cell = cell(length(neighbors_i),1);

%     for j = neighbors_i'
%       Acell{j,1} = param.A_Ni{j};
%     end
%     maxLength = max(cellfun(@length, Acell))
%     Acell_padded = cellfun(@(x) [x zeros(param.ni,maxLength-size(x,2))], ...
%                                  Acell, 'un', 0)
%     A = cell2mat(Acell_padded);
    constraints = [constraints, zDec(setdiff(1:param.nb_subsystems, neighbors_i))...
                                == 0];
    for k = neighbors_i'
        for j = neighbors_i'
           if j <= length(param.Wij{k}) && ~isempty(param.Wij{k}{j})
            Aij{k, j} = param.A_Ni{k}*param.Wij{k}{j}' ;
           else
               Aij{k,j} = zeros(param.ni);
           end 
        end
        Plift_cell{k} = zDec(k)*param.Pi{k};
    end
    A = cell2mat(Aij);
    %Acl = (A+Bblk*Kblk);
    Acl = (A+BiLift*param.K_Ni{i});
    PiLift = blkdiag(Plift_cell{:});
    idx_i_mat = diag(neighbors_i == i);
    QiLift = kron(idx_i_mat, Qi);
    LMI_i = PiLift - QiLift - param.K_Ni{i}'*Ri*param.K_Ni{i} - ...
            Acl'*Pblk*Acl;
%     shiftDiag = 0;
%     for k=1:size(LMI_i,1)
%         constraints = [constraints, LMI_i(k, i*param.ni-1+shiftDiag) >= ...
%                       sum(abs(LMI_i(k,:)))-abs(LMI_i(k, i*param.ni-1+shiftDiag))];
%         shiftDiag = shiftDiag+1;
%     end
    constraints = [constraints, Qi >= 1e-5*eye(param.ni), ...
                                Ri >= 1e-5 * eye(param.nu)];
    constraints = [constraints, LMI_i>=0];
%    
    objective = sum(zDec(neighbors_i).^2);
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0); %options
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
        sprintf('MOSEK solver thinks it is infeasible for system %d', i)
    end
    
    Qi = value(Qi);
    Ri = value(Ri);
    decVariables = value(zDec);
end