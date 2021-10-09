%% main script
clear
% USE MOSEK as solve
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'
param = base_parameters;
K = controller_passivity(param.A_1, param.B_1, param.C_1, param.F_1, ...
                         param.U, param.W);