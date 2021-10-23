%% Algorithm 1: Offline distributed MPC synthesis
% Author:
%   Nicolas Hoischen
% BRIEF:

function [X,U] = mpc_online(x0, alpha, Q_Ni, Ri, Pi, Gamma_i)
param = param_coupled_oscillator;
length_simulation = 10;
M = param.number_subsystem;
x = x0;
for n = 1:length_simulation
    [u0, x_Ni] = first_input(x, alpha, Q_Ni, Ri, Pi);
    x = param.global_sysd.A *x + param-global_sysd.B*u0;
    for i=M:-1:1
        alpha(i) = alpha(i) + x_Ni{i}'*Gamma_i{i}*x_Ni{i};
    end
end

function [U0,x_Ni]  = first_input(x0, alpha_i, Q_Ni, Ri, Pi)

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
X = sdpvar(repmat(nx,1,N),ones(1,N),'full');

% CONSTRAINTS FOR EACH TIME STEP, GLOBAL CENTRALIZED SYSTEM
for k = 1:N-1
    % system dynamics constraints
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
                    + (Z{i}*u{N})'*Ri{i}*(Z{i}*u{N}); % uf = Z{i}*U{i}
    % sum of Vf
    objective = objective + (param.U{i}*X{N})'*Pi{i}*(param.U{i}*X{N}); 
    constraints = [constraints, (param.U{i}*X{N})'*Pi{i}*(param.U{i}*X{N})...
                                <= alpha(i)];
end
% parameter for initial condition
constraints = [constraints, X{1} == x0];
% generate yalmip optimizer object
ops = sdpsettings('verbose',0,'solver','quadprog'); %options
U0 = optimizer(constraints,objective,ops,x0,u{1});
for i = M:-1:1
    x_Ni{i} = param.W{i}*value(X{N});
end
end