%% Main script
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   Interface to test the passive controller
%%
clear
% USE MOSEK as solver (ADD to path)
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'

%% Test coupled oscillator
param = param_coupled_oscillator; % load parameters
%param = param_2_DGU;
% Non-distributed system parameters
Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
Cd = param.global_sysd.C;
x0 = [0.3; 0.5; 0.2; 0.3];
length_sim = 30;
control_type = "Passivity";
[X,U] = simulate_system(@controller_passivity, x0,length_sim, control_type, param);
config = "GENERAL";
sprintf("cost with passive controller is equal to %d", compute_QR_cost(X,U,Q,R))
plot_states_coupled_osci(X,U,config, control_type);
control_type = "LQR";
Q = 100*eye(size(Ad)); R = eye(size(Bd,2));
[X,U, Pinf] = simulate_system(@controller_passivity, x0,length_sim, control_type,...
                        param,Q,R);
sprintf("cost with lqr controller is equal to %d", x0'*Pinf*x0)
sprintf("%d", compute_QR_cost(X,U,Q,R))
plot_states_coupled_osci(X,U,config, control_type);

% elseif param.name == "2_DGU" 
%     sysd_pass = ss(Ad+Bd*Kpassive, Bd,Cd, [], param.Ts);
%     % Comparison to LQR Control
%     Q = 100*eye(6); R = eye(2);
%     Klqr= dlqr(Ad, Bd, Q,  R);
%     sysd_lqr = ss(Ad - Bd * Klqr,[],Cd, [],param.Ts);
% 
%     t = 0:param.Ts:10;
%     u = [-param.Vr1/param.Vin_1+ K1*[-param.Vr1 0 -param.Vr1]'; ...
%          -param.Vr2/param.Vin_2+ K2*[-param.Vr2 0 -param.Vr2]'];
%     u = repmat(u,1,size(t,2))' ;
%     lsim(sysd_pass,u,t) % u matrix with dimensions Nt by Nu
%  end   
%  
 
 %% Test 2; Offline Distributed synthesis coupled oscillator
 clear
 param = param_coupled_oscillator;
 Q_Ni = {}; Ri = {};
 for i = 1:param.number_subsystem
     m_Ni = size(param.W{i},1);
     Q_Ni{i} =100*eye(m_Ni);
     Ri{i} = 1*eye(size(param.Bi{i},2));
 end
 use_passivity = true; 
 [P, Gamma_Ni, alpha_i] = offline_distributed_MPC(Q_Ni, Ri, "COUPLED_OSCI", ...
                                                  use_passivity);
alpha = zeros(param.number_subsystem,1);
 for i=1:param.number_subsystem
    alpha(i) = alpha_i; % same alpha for every subsystem in the beginning
 end
 % send to online controller
 x0 = [0.1 -0.1; 0 -0]; % columns are subsystem i
length_sim = 50;
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
    plot(position{1}, 'r-+');
    plot(position{2}, 'b-*');
    legend("mass 1", "mass 2");
    grid on
    hold off
    subplot(2,1,2)
    title('velocities');
    hold on
    plot(velocity{1}, 'r-+');
    plot(velocity{2}, 'b-*');
    legend("mass 1", "mass 2");
    grid on
    hold off

    figure()
    title("Controller  " + control_type)
    hold on
    plot(controller(1,:), 'r-+');
    plot(controller(2,:), 'b-*');
    legend("U1", "U2");
    grid on
    hold off
end
 
function cost = compute_QR_cost(X,U,Q,R)
cost = 0;
    for i=1:length(X)-1
        cost = cost + X{i}'*Q*X{i} + U{i}'*R*U{i};
    end
    
end