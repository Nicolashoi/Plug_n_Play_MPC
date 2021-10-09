function Ki = controller_passivity(A, B, C, F, U, W)
    E = sdpvar(2,2, 'full');
    H = diag(sdpvar(2,1)); 
    K = sdpvar(1,2, 'full');
    G = K*E;
    S = diag(sdpvar(1));
    % equation 7
    LMI = [E, 1/2*E*C', (A*E + B*G)', E;...
           1/2*C*E, 1/2*S + 1/2*S', F', zeros(size(F',1),size(E,2));...
           (A*E+B*G), F, E, zeros(size(F,1), size(E,2));...
           E, zeros(size(E,1), size(F,2)), zeros(size(E,1), size(E,2)), H];
    constraints = [];
    constraints = [LMI >= 0, E >= eye(2), H >= 0];
    
    %Theorem 1
    epsilon_i = sdpvar(1);
    constraints = [constraints, epsilon_i >= 0.1, E >= epsilon_i*eye(2)];
    epsilon_0 = 0.5;
    for j=1:size(H,1)
        constraints = [constraints, H(j,j) <= 1/(norm(W(j,:),1)+ epsilon_0)];
    end
    for k=1:size(S,1)
       constraints = [constraints, S(k,k) <= 1/norm(U(k,:))]; 
    end
    diagnostics = optimize(constraints, [], sdpsettings('solver', 'MOSEK'))
    if diagnostics.problem == 0
     disp('Solver thinks it is feasible')
    elseif diagnostics.problem == 1
     disp('Solver thinks it is infeasible')
    else
     disp('Something else happened')
    end
    Ki = value(K);
end