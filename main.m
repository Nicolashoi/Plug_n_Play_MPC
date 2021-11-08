%% Main script
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   Interface to test the passive controller
%%
clear
% USE MOSEK as solver (ADD to path)
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'

%% Choose SYSTEM TO USE
% 1: coupled oscillator with spring only, no dampers
% 2: coupled damped oscillator, spring + oscillators
% 3: DGU units not well implemented yet
clear
close
system = 3;

%% TEST 1: PASSIVITY VS LQR CONTROLLER
close all
param = choose_system(system);
% Non-distributed system parameters
Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
Cd = param.global_sysd.C;

length_sim = 500;

% USING PASSIVITY (DISTRIBUTED CONTROL)
control_type = "Passivity";
config = "DISTRIBUTED";
x0_1 = [0.2; 0.1; 0]; x0_2 = [-0.2; 0.1; 0];
x0 = [x0_1 x0_2]; % columns are subsystem i
[X,U] = simulate_system(@controller_passivity, x0,length_sim, control_type, param);

Q = 100*eye(size(Ad)); R = eye(size(Bd,2));
sprintf("cost with passive controller is equal to %d", compute_QR_cost(X,U,Q,R,config))
if param.name == "COUPLED_OSCI"
    plot_states_coupled_osci(X,U,config, control_type);
end

% USING LQR CONTROL
control_type = "LQR";
config = "GENERAL";
x0 = [0.2; 0.1; -0.2; 0.1];
[X,U, Pinf] = simulate_system(@controller_passivity, x0,length_sim, control_type,...
                        param,Q,R);
sprintf("cost of lqr controller using closed form solution %d", x0'*Pinf*x0)
sprintf("cost of lqr finite using function %d", compute_QR_cost(X,U,Q,R, config))
if param.name == "COUPLED_OSCI"
    plot_states_coupled_osci(X,U,config, control_type);
end
 
 %% TEST 2; Offline Distributed synthesis with or without passivity
 % Compute Ki, Pi using passivity or Lyapunov, solve LP to obtain alpha and 
 % run MPC online
param = choose_system(system);
 Q_Ni = {}; Ri = {};
 for i = 1:param.number_subsystem
     m_Ni = size(param.W{i},1);
     Q_Ni{i} =100*eye(m_Ni);
     Ri{i} = 1*eye(size(param.Bi{i},2));
 end
 use_passivity = false; 
 [P, Gamma_Ni, alpha_i] = offline_distributed_MPC(Q_Ni, Ri, param.name, ...
                                                  use_passivity);
alpha = zeros(param.number_subsystem,1);
 for i=1:param.number_subsystem
    alpha(i) = alpha_i; % same alpha for every subsystem in the beginning
 end
 % send to online controller
 x0 = [0.3 -0.4; 0 -0]; % columns are subsystem i
length_sim = 500;
control_type = "MPC";
[X,U] = simulate_system(@mpc_online, x0,length_sim, control_type, param,...
                         Q_Ni, Ri, P, Gamma_Ni, alpha);
config = "DISTRIBUTED";
plot_states_coupled_osci(X,U, config, control_type)

%% FUNCTIONS 
function plot_states_coupled_osci(X,U, config, control_type)  
    if config == "GENERAL"
        states = cell2mat(X);
        position{1} = states(1:4:end); position{2} = states(3:4:end);
        velocity{1} = states(2:4:end); velocity{2} = states(4:4:end); 
        controller = cell2mat(U');%first row u1, second row u2, column are timesteps
    elseif config == "DISTRIBUTED"
        states = cell2mat(X);
        position{1} = states(1:2:end,1);
        position{2} = states(1:2:end,2);
        velocity{1} = states(2:2:end,1);
        velocity{2} = states(2:2:end,2);
        controller = cell2mat(U)'; %first row u1, second row u2, column are timesteps
    else
        error("not implemented configuration in plot states");
    end
        
    figure()
    sgtitle(control_type);
    subplot(2,1,1)
    title('Positions');
    hold on
    plot(position{1}, 'r-');
    plot(position{2}, 'b-');
    legend("mass 1", "mass 2");
    grid on
    hold off
    subplot(2,1,2)
    title('velocities');
    hold on
    plot(velocity{1}, 'r-');
    plot(velocity{2}, 'b-');
    legend("mass 1", "mass 2");
    grid on
    hold off

    figure()
    title("Controller  " + control_type)
    hold on
    plot(controller(1,:), 'r-');
    plot(controller(2,:), 'b-');
    legend("U1", "U2");
    grid on
    hold off
end
 
function param = choose_system(system)
    switch system
        case 1
            param = param_coupled_oscillator;
        case 2
            param = param_mass_damper;
        case 3
            param = param_2_DGU;
        otherwise
            error("system not implemented yet")
    end

end

function cost = compute_QR_cost(X,U,Q,R, config)
    cost = 0;
    if config == "GENERAL"
        for i=1:length(X)-1
            cost = cost + X{i}'*Q*X{i} + U{i}'*R*U{i};
        end
    elseif config == "DISTRIBUTED"
        for i=1:length(X)-1
            x = reshape(X{i},[],1); u = U{i}';
            cost = cost + x'*Q*x + u'*R*u;
        end
    else
        error("wrong configuration, choose between general or distributed");
    end   
end