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
    A_2 = [0, 1/Cap2, 0; -1/L2, -R2/L2, 0; 1/Ts, 0, 0];
    B_2 = [0; Vin_2/L2; 0];
    F_2 = [1/Cap2; 0; 0];
    
    %% Laplacian
    l11 = 0; l12 = 1/R12; l21 = l12; l22 = 0;
    L = [l11+l12,-l12; -l21; l21 + l22];
    


end