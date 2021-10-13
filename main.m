%% Main script
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   Interface to test the passive controller
%%
clear
% USE MOSEK as solver
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'

%% Test coupled oscillator
param = param_coupled_oscillator; % load parameters
[K1, D1, P1, Gamma_1] = controller_passivity(param.A_1, param.B_1,...
                                             param.C_1, param.F_1, ...
                                             param.U, param.W);
[K2,D2, P2, Gamma_2] = controller_passivity(param.A_2, param.B_2,...
                                            param.C_2, param.F_2, ...
                                            param.U, param.W);
% Construct subsystems
sys1d = ss(param.A_1+param.B_1*K1, param.F_1, param.C_1, D1, -1);
sys2d = ss(param.A_2+param.B_2*K2, param.F_2, param.C_2, D2, -1);
disp("Is subsystem 1 passive ?"); disp(isPassive(sys1d));
disp("Is subsystem 2 passive ?"); disp(isPassive(sys2d));  

% Non-distributed system parameters
Ad = param.global_sysd.A; Bd = param.global_sysd.B; Cd = param.global_sysd.C;
Gamma = blkdiag(Gamma_1, Gamma_2); D = blkdiag(D1, D2); % definitions
epsilon_0 = 0.5; % same as in controller passivity.m
LMI_15 = [Gamma - epsilon_0*eye(size(Gamma)) + Cd'*param.L_tild'*Cd, Cd'*param.L_tild';...
          param.L_tild'*Cd, inv((D+D')/2)];
if all(eig(LMI_15)) > 0
    disp("Eq.15 verified, global system is A.S according to Lemma 2");
end
disp("Minimum eigenvalue of the dissipation rate matrix Gamma");
disp(min(eig(Gamma)));
% Construct global system with the subsystems controllers
Kpassive = [K1, zeros(size(K1,1), size(K2,2)); zeros(size(K2,1), size(K1,2)), K2];
sysd_pass = ss(Ad+Bd*Kpassive, [],Cd, [], -1);

% Comparison to LQR Control
 Q = 100*eye(4); R = eye(2);
 Klqr = dlqr(Ad, Bd, Q,  R);
 sysd_lqr = ss(Ad - Bd * Klqr,[],Cd, [],-1);
 
 % Plot both controllers for initial conditions
 x0 = [1; 0.5; 0; 0.3];
 initial(sysd_lqr, x0)
 hold on
 initial(sysd_pass, x0)
 legend("LQR Control", "Passive Controller");
 grid on
 hold off
 
 %% Test 2 distributed connected DGU's
 param = param_2_DGU;
 [K1, D1, P1, Gamma_1] = controller_passivity(param.A_1, param.B_1,...
                                             param.C_1, param.F_1, ...
                                             param.U, param.W);
 [K2,D2, P2, Gamma_2] = controller_passivity(param.A_2, param.B_2,...
                                            param.C_2, param.F_2, ...
                                            param.U, param.W);
% Construct subsystems
sys1d = ss(param.A_1+param.B_1*K1, param.F_1, param.C_1, D1, -1);
sys2d = ss(param.A_2+param.B_2*K2, param.F_2, param.C_2, D2, -1);
disp("Is subsystem 1 passive ?"); disp(isPassive(sys1d));
disp("Is subsystem 2 passive ?"); disp(isPassive(sys2d));  



 