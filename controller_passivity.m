%% YALMIP & MOSEK Semi-definite program
% Author:
%   Nicolas Hoischen
% BRIEF:
%   Controller function to compute Ki to make the subsystem passive. 
%   From "Passivity-Based Decentralized Control for Discrete-Time Large-Scale
%         Systems"
% INPUT:
%   Ai, Bi, Ci, Fi: LTI dynamics of subsystem i
%   U: L_tilde*C
%   W: C^T * L_tilde where L_tilde is the augmented laplacian
% OUTPUT:
%   Ki: decentralized passivity-based control state feedback gain
%   Di: positive-definite diagonal matrix
%   Pi: positive-definite matrix
%   Gamma_i: positive-definite diagonal matrix

function [Ki, Di, Pi, Gamma_i] = controller_passivity(A, B, C, F, L_tilde,...
                                                      C_global)
    %% Interconnection Variables
    U = L_tilde*C_global;
    W = C_global'*L_tilde;

    ni = size(A,1);
    mi = size(B,2);
    constraints = [];
    objective = 0;%trace(H);
    
    %% Decision Variables
    E = sdpvar(ni, ni); %symmetric P.D
    H = diag(sdpvar(ni,1)); % diagonal matrix as defined in Th.1
    G = sdpvar(mi,ni, 'full');
    S = diag(sdpvar(mi,1)); % diagonal matrix as defined in Th.1
   
    %% Equation 7
    LMI = [E, 1/2*E*C', (A*E + B*G)', E;...
           1/2*C*E, 1/2*S + 1/2*S', F', zeros(size(F',1),size(E,2));...
           (A*E+B*G), F, E, zeros(size(F,1), size(E,2));...
           E, zeros(size(E,1), size(F,2)), zeros(size(E)), H];
    % add to constraints   
    constraints = [constraints, LMI >= 0]; %0-1e-2*eye(size(LMI,1))]; 

    %% Theorem 1
    %epsilon_i = sdpvar(1); %define value
    epsilon_i = 1e-5;
    constraints = [constraints, E >= epsilon_i*eye(ni), ...
                   H >= epsilon_i*eye(ni), S >= epsilon_i*eye(mi)];
    epsilon_0 = 1e-5;
    for j=1:size(H,1)
        constraints = [constraints, H(j,j) <= 1/(norm(W(j,:),1)+ epsilon_0)];
    end
    for k=1:size(S,1)
       if norm(U(k,:),1) ~= 0
            constraints = [constraints, S(k,k) <= 1/norm(U(k,:),1)]; 
       end
    end
    %% OPTIMIZER
    ops = sdpsettings('solver', 'MOSEK', 'verbose',0);
    %ops.mosek.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON';
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem == 1
        disp('MOSEK solver thinks it is infeasible')
    end
    %% MAP
    Pi = inv(value(E)); Ki = value(G)*Pi;
    Gamma_i = inv(value(H)); Di = value(S);
end