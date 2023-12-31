%% Semi-definite program to compute passivating feedback controllers
% Author:
%   Nicolas Hoischen
% BRIEF:
%   Controller function to compute Ki to make each subsystem strictly passive. 
%   Intended for the redesign phase based on passivity.
%   From Ahmed Aboudonia et al. "Passivity-Based Decentralized Control for 
%   "Discrete-Time Large-Scale Systems"
%   IEEE Control Systems Letters 5.6 (2021)
% INPUT:
    % Ai, Bi, Ci, Fi: LTI dynamics of subsystem i
    % L_tilde: Graph augmented Laplacian
    % C_global: C matrix of the overall system
% OUTPUT:
    % Ki: decentralized passivity-based control state feedback gain
    % Di: positive-definite diagonal matrix
    % Pi: positive-definite matrix
    % Gamma_i: positive-definite diagonal matrix
    % feasible: boolean true if problem is feasible

%%
function [Ki, Di, Pi, Gamma_i, feasible] = controller_passivity(A, B, C, F, L_tilde,...
                                                      C_global, i)
    %% Interconnection Variables
    U = L_tilde*C_global;
    W = C_global'*L_tilde';
    ni = size(A,1);
    mi = size(B,2);
    Ui = U(i,:);
    Wi = W(i:i+ni-1,:);
    constraints = [];
     
    %% Decision Variables
    E = sdpvar(ni, ni); %symmetric P.D
    H = diag(sdpvar(ni,1)); % diagonal matrix as defined in Th.1
    G = sdpvar(mi,ni, 'full');
    S = diag(sdpvar(mi,1)); % diagonal matrix as defined in Th.1
    objective = trace(H); %norm(G);%trace(H)-trace(E);% yalmip always assumes minimization
    %% Equation 7
    LMI = [E, 1/2*E*C', (A*E + B*G)', E;...
           1/2*C*E, 1/2*S + 1/2*S', F', zeros(size(F',1),size(E,2));...
           (A*E+B*G), F, E, zeros(size(F,1), size(E,2));...
           E, zeros(size(E,1), size(F,2)), zeros(size(E)), H];
    % add to constraints   
    constraints = [constraints, LMI >= 0]; 

    %% Theorem 1
    %epsilon_i = sdpvar(1); %define value
    epsilon_i = 1e-5;
    constraints = [constraints, E >= epsilon_i*eye(ni), ...
                   H >= epsilon_i*eye(ni), S >= epsilon_i*eye(mi)];
    epsilon_0 = 1e-5;
    for j=1:ni
       constraints = [constraints, H(j,j) <= 1/(norm(Wi(j,:),1)+ epsilon_0)];
    end
    for k=1:mi
       if norm(U(k,:),1) ~= 0
           constraints = [constraints, S(k,k) <= 1/norm(Ui(k,:),1)]; 
       else 
           disp("Warning: division by 0 in controller passivity Theorem 1");
       end
    end
    %% OPTIMIZER
    feasible = 1;
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0);
    %ops.mosek.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON';
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
       sprintf('MOSEK solver thinks it is infeasible for system %d', i)
       feasible = 0;
    end
    %% MAP
    Pi = inv(value(E)); Ki = value(G)*Pi;
    Gamma_i = inv(value(H)); Di = value(S);
end