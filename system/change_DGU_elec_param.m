%% Change DGU parameters to create a *mat file
clear
% System 1
R(1) = 3e-3; Li(1) = 92e-6; Ct(1) = 75e-6; 
Vin(1) = 100; 

% System 2
R(2) = 1.5e-3; Li(2) = 102e-6; Ct(2) = 50e-6;
Vin(2) = 100;

% System 3
R(3) = 1.7e-3; Li(3) = 100e-6; Ct(3) = 56e-6;
Vin(3) = 100;

% System 4
R(4) = 1.6e-3; Li(4) = 94e-6; Ct(4) = 62e-6;
Vin(4) = 100;

% System 5
R(5) = 1.5e-3; Li(5) = 92e-6; Ct(5) = 69e-6;
Vin(5) = 100;

% System 6
R(6) = 1.3e-3; Li(6) = 94e-6; Ct(6) = 62e-6;
Vin(6) = 100;

% Rij link
R12 = 1.75; R23 = 3.5; R24 = 1.75; R35 = 2.8; R36 = 1.75;

Rij = [0 R12 0 0 0 0; R12 0 R23 R24 0 0; 0 R23 0 0 R35 R36; ...
       0 R24 0 0 0 0; 0 0 R35 0 0 0; 0 0 R36 0 0 0];
Agraph = [0 1/R12 0 0 0 0; 1/R12 0 1/R23 1/R24 0 0; 0 1/R23 0 0 1/R35 1/R36; ...
       0 1/R24 0 0 0 0; 0 0 1/R35 0 0 0; 0 0 1/R36 0 0 0];
nb_subsystems = 6; 
size_subsystem = 2;
L = diag(sum(Agraph))-Agraph;
% save to mat file
save('DGU_electrical_param.mat');