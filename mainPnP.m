%% Main File for PnP
clear
% USE MOSEK as solver (ADD to path)
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'
addpath(genpath(cd));
%% Initial DGU Network
clear
utils = utilityFunctions;
PnP = operationPnP;
filename = 'config_DGU_1.txt';
[nb_subsystems, Vin,R,L,C, Vmax, Vmin, Imax, Imin] = utils.importData(filename);

Rij_mat = zeros(nb_subsystems);
Rij_mat(1,2) = 1.75; Rij_mat(2,3) = 3.5; Rij_mat(2,4) = 1.75; 
Rij_mat(3,5) = 2.8; %Rij_mat(3,6) = 2; 
Rij_mat = Rij_mat + tril(Rij_mat',1);
%% Instantiate DGU NETWORK
dguNet = DGU_network(nb_subsystems);
Vr = linspace(49.95, 50.2, nb_subsystems);% references
Il = linspace(2.5, 7.5, nb_subsystems);
% set Electrical parameters and Dynamics for ALL the subsystems in the network
for i=1:nb_subsystems
    dguNet = dguNet.initElecParam(i,Vin(i), Vr(i), Il(i), R(i), C(i), L(i), ...
                                  Vmax(i), Vmin(i), Imax(i), Imin(i));
end

%% First SCENARIO
activeDGU_scen1 = 1:1:5; % Initially 5 DGU are active
dguNet = dguNet.setActiveDGU(Rij_mat, activeDGU_scen1);
figure()
plot(dguNet.NetGraph, 'EdgeLabel', dguNet.NetGraph.Edges.Weight, 'Marker', 's', 'NodeColor','r', ...
      'MarkerSize', 7);
title("Initial Scenario: 5 DGUs active");
dguNet = dguNet.initDynamics();
delta_config = false;
dguNet = dguNet.setConstraints(delta_config);
config = "DISTRIBUTED";
control_type = "MPC online";
[x0, Q_Ni, Ri] = utils.tuningParam(dguNet, delta_config);
length_sim = 40;
dguNet = PnP.setPassiveControllers(dguNet);
% Test reconf. terminal ingredients MPC with initial DGU network
[X, U] = PnP.mpc_DGU_tracking(@mpc_online_2, x0, length_sim, dguNet, Q_Ni, Ri);
dguNet.plot_DGU_system(X,U, config, control_type, dguNet); % plot results

%% Scenario 2: Add a DGU to the network: connect DGU 6 to DGU 3
dguPos = 6;
activeDGU_scen2 = 1:1:6;
Rij_mat(3,6) = 0.75; Rij_mat(6,3) = Rij_mat(3,6);  % New link
dguNet2 = dguNet; % copy instance class to new object (keep both objects to compare)
dguNet2 = dguNet2.setActiveDGU(Rij_mat, activeDGU_scen2);
dguNet2 = dguNet2.initDynamics(); % recompute Dynamics (changed with integration of DGU 6)
plot(dguNet2.NetGraph, 'EdgeLabel', dguNet2.NetGraph.Edges.Weight, 'Marker', 's', 'NodeColor','r', ...
      'MarkerSize', 7);
dguNet2 = dguNet2.setConstraints(delta_config);

dguNet2 = PnP.redesignPhase(dguNet2, dguNet2.NetGraph,dguPos, "add");
[~, Q_Ni, Ri] = utils.tuningParam(dguNet2, delta_config);
for i=activeDGU_scen1
    x0_scen2{i} = X{end}(:,i); % previous states start
end
x0_scen2{dguPos} = [50;0];
[Xscen2,Uscen2] = PnP.mpc_DGU_tracking(@mpc_online_2, x0_scen2, length_sim, dguNet2, Q_Ni, Ri);
dguNet2.plot_DGU_system(Xscen2, Uscen2, config, control_type, dguNet2); % plot results

%% Remove DGU 4
dguDelete = 4;
Rij_mat(dguDelete,:) = 0; Rij_mat(:,dguDelete) = 0;
activeDGU_scen3 = [1 2 3 5 6]; 
dguNet3 = dguNet2;
dguNet3 = dguNet3.setActiveDGU(Rij_mat, activeDGU_scen3);
dguNet3 = dguNet3.initDynamics();
plot(dguNet3.NetGraph, 'EdgeLabel', dguNet3.NetGraph.Edges.Weight, 'Marker', 's', 'NodeColor','r', ...
      'MarkerSize', 7);
dguNet3 = dguNet3.setConstraints(delta_config);

dguNet3 = PnP.redesignPhase(dguNet3, dguNet2.NetGraph, dguDelete, "delete");
for i=activeDGU_scen2
    x0_scen3{i} = Xscen2{end}(:,i); % previous states start
end
%[~, Q_Ni, Ri] = utils.tuningParam(dguNet3, delta_config);
[Xscen3,Uscen3] = PnP.mpc_DGU_tracking(@mpc_online_2, x0_scen3, length_sim, dguNet3, Q_Ni, Ri);
dguNet3.plot_DGU_system(Xscen3, Uscen3, config, control_type, dguNet3); % plot results