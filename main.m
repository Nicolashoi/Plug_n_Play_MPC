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

% Compute a passive controller for each subsystem
[K1, D1, P1, Gamma_1] = controller_passivity(param.A_1, param.B_1,...
                                             param.C_1, param.F_1, ...
                                             param.L_tilde, Cd);
[K2,D2, P2, Gamma_2] = controller_passivity(param.A_2, param.B_2,...
                                            param.C_2, param.F_2, ...
                                            param.L_tilde,Cd);
% Controlled subsystems
sys1d = ss(param.A_1+param.B_1*K1, param.F_1, param.C_1, D1, -1);
sys2d = ss(param.A_2+param.B_2*K2, param.F_2, param.C_2, D2, -1);
disp("Is subsystem 1 passive ?"); disp(isPassive(sys1d));
disp("Is subsystem 2 passive ?"); disp(isPassive(sys2d));  

% Check if conditions for A.S are satisfied
Gamma = blkdiag(Gamma_1, Gamma_2); D = blkdiag(D1, D2); % definitions
epsilon_0 = 1e-5; % same as in controller passivity.m
LMI_15 = [Gamma - epsilon_0*eye(size(Gamma)) + Cd'*param.L_tilde'*Cd, Cd'*param.L_tilde';...
          param.L_tilde'*Cd, inv((D+D')/2)];
if all(eig(LMI_15)) > 0
    disp("Eq.15 verified, global system is A.S according to Lemma 2");
end
disp("Minimum eigenvalue of the dissipation rate matrix Gamma");
disp(min(eig(Gamma)));

% Construct global system with the subsystems controllers
Kpassive = blkdiag(K1,K2);
if param.name == "COUPLED_OSCI"
    sysd_pass = ss(Ad+Bd*Kpassive, [],Cd, [], param.Ts);
    Q = 100*eye(4); R = eye(2);
    Klqr= dlqr(Ad, Bd, Q,  R);
    sysd_lqr = ss(Ad - Bd * Klqr,[],Cd, [],param.Ts);
    x0 = [1; 0.5; 0; 0.3];
    initial(sysd_lqr, x0) 
    % Plot both controllers for initial conditions
    hold on
    initial(sysd_pass, x0)
    legend("LQR Control", "Passive Controller");
    grid on
    hold off
elseif param.name == "2_DGU" 
    sysd_pass = ss(Ad+Bd*Kpassive, Bd,Cd, [], param.Ts);
    % Comparison to LQR Control
    Q = 100*eye(6); R = eye(2);
    Klqr= dlqr(Ad, Bd, Q,  R);
    sysd_lqr = ss(Ad - Bd * Klqr,[],Cd, [],param.Ts);

    t = 0:param.Ts:10;
    u = [-param.Vr1/param.Vin_1+ K1*[-param.Vr1 0 0]'; ...
         -param.Vr2/param.Vin_2+ K2*[-param.Vr2 0 0]'];
    u = repmat(u,1,size(t,2))' ;
    lsim(sysd_pass,u,t) % u matrix with dimensions Nt by Nu
 end   
 
 
 %% Test 2; Offline Distributed synthesis coupled oscillator
 param = param_coupled_oscillator;
 Q_Ni = {}; Ri = {};
 for i = 1:param.number_subsystem
     m_Ni = size(param.W{i},1);
     Q_Ni{i} =100*eye(m_Ni);
     Ri{i} = 1*eye(size(param.B_1,2));
 end
 [P, Gamma_Ni, alpha_i] = offline_distributed_MPC(Q_Ni, Ri, "COUPLED_OSCI");



 