%% Main script
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   Interface to test the passive controller
%%
clear
% USE MOSEK as solver (ADD to path)
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'
addpath(genpath(cd));
%% Choose SYSTEM TO USE
% 1: coupled oscillator with spring only, no dampers
% 2: coupled damped oscillator, spring + oscillators
% 3: DGU units not well implemented yet
% 4: DGU units in delta form
clear
close
system = 4;

%% TEST 1: PASSIVITY VS LQR CONTROLLER
close all
[param, x0] = choose_system(system);
% Non-distributed system parameters
Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
Cd = param.global_sysd.C;

length_sim = 50;
utils = utilityFunctions;
% USING PASSIVITY (DISTRIBUTED CONTROL)
control_type = "PASSIVITY";
config = "DISTRIBUTED";
[X,U] = simulate_system(@controller_passivity, x0,length_sim, control_type, param);

Q = eye(size(Ad)); R = eye(size(Bd,2));
sprintf("cost with passive controller is equal to %d", ...
         utils.compute_QR_cost(X,U,Q,R,config))
plot_simulation(X,U, config, control_type, param, utils)

% USING LQR CONTROL
control_type = "LQR";
config = "GENERAL";

[X,U, Pinf] = simulate_system(@controller_passivity, x0,length_sim, control_type,...
                        param,Q,R);
LQR_cost = vertcat(x0{:})'*Pinf*vertcat(x0{:});
sprintf("cost of lqr controller using closed form solution %d", LQR_cost)
sprintf("cost of lqr finite using function %d",...
        utils.compute_QR_cost(X,U,Q,R, config))
plot_simulation(X,U, config, control_type, param, utils)
 
 %% TEST 2; Offline Distributed synthesis with or without passivity
 % Compute Ki, Pi using passivity or Lyapunov, solve LP to obtain alpha and 
 % run MPC online
[param, x0] = choose_system(system);
utils = utilityFunctions;
Q_Ni = {}; Ri = {};
 for i = 1:param.number_subsystem
     m_Ni = size(param.W{i},1);
     Q_Ni{i} =1*eye(m_Ni);
     Ri{i} = 1*eye(size(param.Bi{i},2));
 end
 use_passivity = false; 
 [P, Gamma_Ni, alpha_i] = offline_distributed_MPC(Q_Ni, Ri, param, ...
                                                  use_passivity);
alpha = zeros(param.number_subsystem,1);
 for i=1:param.number_subsystem
    alpha(i) = alpha_i; % same alpha for every subsystem in the beginning
 end
 % send to online controller
length_sim = 30;
control_type = "MPC";
[X,U] = simulate_system(@mpc_online, x0,length_sim, control_type, param,...
                         Q_Ni, Ri, P, Gamma_Ni, alpha);
config = "DISTRIBUTED";
plot_simulation(X,U, config, control_type, param, utils)

%% TEST 3
[param, x0] = choose_system(system);
utils = utilityFunctions;
Q_Ni = {}; Ri = {};
 for i = 1:param.number_subsystem
     m_Ni = size(param.W{i},1);
     Q_Ni{i} =1*eye(m_Ni);
     Ri{i} = 1*eye(size(param.Bi{i},2));
 end
 % send to online controller
length_sim = 30;
control_type = "MPC";
[X,U] = simulate_system(@mpc_online_2, x0,length_sim, control_type, param,...
                         Q_Ni, Ri, Pi, Gamma_Ni, alpha);
config = "DISTRIBUTED";
plot_simulation(X,U, config, control_type, param, utils)
%% FUNCTIONS 

function plot_simulation(X,U, config, control_type, param, utils)
    if param.name == "COUPLED_OSCI"
        utils.plot_states_coupled_osci(X,U,config, control_type, param);
    elseif param.name == "2_DGU" || param.name == "DGU_delta" 
        utils.plot_DGU_system(X,U, config, control_type, param);
    end
end


function [param, x0] = choose_system(system)
    switch system
        case 1
            param = param_coupled_oscillator;
            x0{1} = [0.3; 0.1]; 
            x0{2} = [-0.2; 0.1];
        case 2
            param = param_mass_damper;
            x0{1} = [0.3; 0.1]; 
            x0{2} = [-0.2; 0.1];
        case 3
            param = param_2_DGU;
            x0{1} = [50.1; 0; 0]; 
            x0{2} = [49.9; 0 ;0];
        case 4
            param = param_DGU_delta;
            x0 = cell(1, param.number_subsystem);
            for i=1:param.number_subsystem
                x0{i} = [50;5];
            end
        otherwise
            error("system not implemented yet")
    end

end