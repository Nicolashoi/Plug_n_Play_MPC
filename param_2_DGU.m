%% Parameter file for a Network of 2 DGUs %%
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   2 DGU connected by a resistor

function param = param_2_DGU
    %% Converter Parameters
    R1 = 0.2; L1 = 1.8e-3; Cap1 = 2.2e-3;
    R2 = 0.3; L2 = 2e-3; Cap2 = 1.9e-3;
    R12 = 0.05; 
    Vin_1 = 48; Vin_2 = 48;
    
    %% System dynamics
    Ts = 1e-5; % sampling time
    A_1 = [0, 1/Cap1, 0; -1/L1, -R1/L1, 0; 1/Ts, 0, 0];
    B_1 = [0; Vin_1/L1; 0];
    F_1 = [1/Cap1; 0; 0];
    C_1 = [1, 0, 0];
    A_2 = [0, 1/Cap2, 0; -1/L2, -R2/L2, 0; 1/Ts, 0, 0];
    B_2 = [0; Vin_2/L2; 0];
    F_2 = [1/Cap2; 0; 0];
    C_2 = C_1;
    %% Laplacian
    l11 = 0; l12 = 1/R12; l21 = l12; l22 = 0;
    L = [l11+l12,-l12; -l21, l21 + l22];
    L_tild = kron(L, eye(size(C_1,1)));
    C = [C_1, zeros(size(C_1,1),size(C_2,2)); zeros(size(C_2,1),size(C_1,2)), C_2];
    U = L_tild*C;
    W = C'*L_tild;
    %% Exact discretization
    sys1 = ss(A_1, [B_1 F_1], [], []);
    sys1_d = c2d(sys1, Ts);
    sys2 = ss(A_2, [B_2 F_2], [], []);
    sys2_d = c2d(sys2, Ts);
    
    %% put everything together
    param.A_1 = sys1_d.A;
    param.B_1 = sys1_d.B(:,1);
    param.F_1 = sys1_d.B(:,2);
    param.A_2 = sys2_d.A;
    param.B_2 = sys2_d.B(:,1);
    param.F_2 = sys2_d.B(:,2);
    param.C_1 = C_1;
    param.C_2 = C_2;
    param.U = U;
    param.W = W;
    param.L_tild = L_tild;

end