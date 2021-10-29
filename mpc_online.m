%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:

function [X,U] = mpc_online(x0, alpha, Q_Ni, Ri, Pi, Gamma_Ni)
param = param_coupled_oscillator;
length_sim = 30;
M = param.number_subsystem;
x = x0;
% X = cell(size(param.global_sysd.A,1),1,length_sim+1);
% U = cell(size(param.global_sysd.B,2),1,length_sim);
X = cell(length_sim+1,1);
U = cell(length_sim,1);
X{1} = x0;
    for n = 1:length_sim
        [u0, x_Ni] = get_input(x, alpha, Q_Ni, Ri, Pi);
        x = param.global_sysd.A *x + param.global_sysd.B*u0;
        X{n+1} = x; U{n} = u0;
        for i=M:-1:1
            alpha(i) = alpha(i) + x_Ni{i}'*Gamma_Ni{i}*x_Ni{i}; %apply the x_Ni we get from system equations after applying x0
        end
    end
end

function [U0,x_Ni]  = get_input(x0, alpha, Q_Ni, Ri, Pi)

    param = param_coupled_oscillator;
    N = 5;
    M = param.number_subsystem;
    A = param.global_sysd.A; B = param.global_sysd.B;
    % create variables for optimizer
    nx = size(A,1);
    nu = size(B,2); 
    objective = 0;
    constraints = [];
    u = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');
    X = sdpvar(repmat(nx,1,N),ones(1,N),'full'); % arrays of nx*1 repeated N times
    X0 = sdpvar(nx,1,'full');
    alpha_var = sdpvar(M,1,'full');
    % CONSTRAINTS FOR EACH TIME STEP, GLOBAL CENTRALIZED SYSTEM
    for k = 1:N-1
        % system dynamics constraints #CHANGE TO DISTRIBUTED SYSTEM A_i,
        % X_i
        constraints = [constraints, X{k+1} == A*X{k} + B*u{k}];
        % state constraints
        constraints = [constraints, param.Gx * X{k+1} <= param.fx];
        % input constraints
        constraints = [constraints, param.Gu * u{k} <= param.fu];
    end
    % FOR EACH SUBSYSTEM
    for i = 1:M
        % sum of local functions l(xi,uf)
        objective = objective + (param.W{i}*X{N})'*Q_Ni{i}*(param.W{i}*X{N})...
                        + (param.Z{i}*u{N-1})'*Ri{i}*(param.Z{i}*u{N-1}); % uf = Z{i}*U{i}
        % sum of Vf
        objective = objective + (param.U{i}*X{N})'*Pi{i}*(param.U{i}*X{N}); 
        constraints = [constraints, (param.U{i}*X{N})'*Pi{i}*(param.U{i}*X{N})...
                                    <= alpha_var(i)];
    end
    % parameter for initial condition
    constraints = [constraints, X{1} == X0];
    ops = sdpsettings('verbose',1); %options
    parameters_in = {X0, alpha_var};
    solutions_out = {[u{:}], [X{:}]};
    MPC_optimizer = optimizer(constraints,objective,ops,parameters_in,solutions_out);
    solutions = MPC_optimizer({x0, alpha});
    U0 = solutions{1}(:,1);
    XN = solutions{2}(:,4);
    for i = M:-1:1
        x_Ni{i} = param.W{i}*XN; % DON?T DO THAT HERE, but in the other function and get x1 not XN
    end
end