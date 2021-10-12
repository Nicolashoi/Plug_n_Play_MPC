%% main script
clear
% USE MOSEK as solver
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'
param = param_coupled_oscillator;
[K1, D1, P1, Gamma_1] = controller_passivity(param.A_1, param.B_1, param.C_1, param.F_1, ...
                         param.U, param.W);
[K2,D2, P2, Gamma_2] = controller_passivity(param.A_2, param.B_2, param.C_2, param.F_2, ...
                         param.U, param.W);
sys1d = ss(param.A_1+param.B_1*K1, param.F_1, param.C_1, D1, -1);
sys2d = ss(param.A_2+param.B_2*K2, param.F_2, param.C_2, D2, -1);
disp("Is subsystem 1 passive ?"); disp(isPassive(sys1d));
disp("Is subsystem 2 passive ?"); disp(isPassive(sys2d));   
Ad = param.global_sysd.A;
Bd = param.global_sysd.B;
Cd = param.global_sysd.C;
Gamma = blkdiag(Gamma_1, Gamma_2);
D = blkdiag(D1, D2);
epsilon_0 = 0.5;
Lemma2 = [Gamma - epsilon_0*eye(size(Gamma)) + Cd'*param.L_tild'*Cd, Cd'*param.L_tild';...
          param.L_tild'*Cd, inv((D+D')/2)];
disp("eigenvalue of 15)");
disp(eig(Lemma2));


%% LQR Controller
 param = param_coupled_oscillator;
 Q = 100*eye(4);
 R = eye(2);
 Ad = param.global_sysd.A;
 Bd = param.global_sysd.B;
 Cd = param.global_sysd.C;
 Klqr = dlqr(Ad, Bd, Q,  R);
 Acl = Ad - Bd * Klqr;
 sysd_lqr = ss(Acl,[],Cd, [],-1);
 Kpassive = [K1, zeros(size(K1,1), size(K2,2)); zeros(size(K2,1), size(K1,2)), K2];
 sysd_pass = ss(Ad+Bd*Kpassive, [],Cd, [], -1);
 initial(sysd_lqr, [1; 0.1;0; 0.3])
 hold on
 initial(sysd_pass, [1; 0.1;0; 0.3])
 hold off
 