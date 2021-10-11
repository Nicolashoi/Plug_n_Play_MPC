%% main script
clear
% USE MOSEK as solve
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'
param = base_parameters;
[K1, D1, P1, Gamma_1] = controller_passivity(param.A_1, param.B_1, param.C_1, param.F_1, ...
                         param.U, param.W);
[K2,D2, P2, Gamma_2] = controller_passivity(param.A_2, param.B_2, param.C_2, param.F_2, ...
                         param.U, param.W);
sys1d = ss(param.A_1+param.B_1*K1, param.F_1, param.C_1, D1, -1);
sys2d = ss(param.A_2+param.B_2*K2, param.F_2, param.C_2, D2, -1);
disp("Is subsystem 1 passive ?"); disp(isPassive(sys1d));
disp("Is subsystem 2 passive ?"); disp(isPassive(sys2d));                  
 %% LQR Controller
 param = base_parameters;
 Q = 100*eye(4);
 R = eye(2);
 Ad = param.global_sysd.A;
 Bd = param.global_sysd.B;
 Cd = param.global_sysd.C;
 Klqr = dlqr(Ad, Bd, Q,  R);
 Acl = Ad - Bd * Klqr;
 sysd = ss(Acl,[],Cd, [],-1);
 initial(sysd, [1; 0.1;0; 0.3])