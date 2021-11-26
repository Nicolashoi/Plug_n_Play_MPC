%% Main File for PnP
clear
% USE MOSEK as solver (ADD to path)
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'
addpath(genpath(cd));
%% Initial DGU Network

filename = 'config_DGU_1.txt';
[nb_subsystems, Vin,R,L,C, Vmax, Vmin, Imax, Imin] = importData(filename);

Rij_mat = zeros(nb_subsystems);
Rij_mat(1,2) = 1.75; Rij_mat(2,3) = 3.5; Rij_mat(2,4) = 1.75; 
Rij_mat(3,5) = 2.8;
Rij_mat = Rij_mat + tril(Rij_mat',1);

dguNet = DGU_network(nb_subsystems, Rij_mat);
% references
Vr = linspace(49.95, 50.2, nb_subsystems);
Il = linspace(2.5, 7.5, nb_subsystems);
% set Electrical parameters
for i=1:nb_subsystems
    dguNet = dguNet.initElecParam(i,Vin(i), Vr(i), Il(i), R(i), C(i), L(i), ...
                                  Vmax(i), Vmin(i), Imax(i), Imin(i));
end
dguNet = dguNet.initDynamics([1:1:nb_subsystems]);
control_type = "MPC_2";
[x0, Q_Ni, Ri, control_type] = tuningParam(dguNet,control_type)
length_sim = 50;
[X,U] = simulate_system(@mpc_online_2, x0,length_sim, control_type, dguNet,...
                         Q_Ni, Ri);
config = "DISTRIBUTED";
dguNet.plot_DGU_system(X,U, config, control_type, dguNet);


function [x0, Q_Ni, Ri, control_type] = tuningParam(dguNet,control_type)
    Q_Ni = cell(1,dguNet.nb_subsystems); 
    Ri = cell(1,dguNet.nb_subsystems);
    x0 = cell(1,dguNet.nb_subsystems);
    for i = 1:dguNet.nb_subsystems
        x0{i} = [50;0]; % second state is Ii - Il
        m_Ni = size(dguNet.W{i},1);
        Q_Ni{i} =1*eye(m_Ni);
        Ri{i} = 1*eye(size(dguNet.Bi{i},2));
    end
end

function [subsystems, Vin, R, L, C, Vmax, Vmin, Imax, Imin] = importData(filename)
    delimiterIn = ';';
    headerlinesIn= 1;
    dataDGU = importdata(filename, delimiterIn, headerlinesIn);
    subsystems = dataDGU.data(1,1);
    Vin = dataDGU.data(:,2);
    R = dataDGU.data(:,3);
    L = dataDGU.data(:,4);
    C = dataDGU.data(:,5);
    Vmax = dataDGU.data(:,6);
    Vmin = dataDGU.data(:,7);
    Imax = dataDGU.data(:,8);
    Imin = dataDGU.data(:,9);
end