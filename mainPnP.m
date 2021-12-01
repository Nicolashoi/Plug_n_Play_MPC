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
Rij_mat(3,5) = 2.8;
Rij_mat = Rij_mat + tril(Rij_mat',1);
%% Instantiate DGU NETWORK
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
%% Test reconf. terminal ingredients MPC with initial DGU network
delta_config = false;
dguNet = dguNet.setConstraints(delta_config);
config = "DISTRIBUTED";
control_type = "MPC online";
[x0, Q_Ni, Ri] = utils.tuningParam(dguNet, delta_config);
length_sim = 40;
dguNet = PnP.setPassiveControllers(dguNet);
[X, U] = PnP.mpc_DGU_tracking(@mpc_online_2, x0, length_sim, dguNet, Q_Ni, Ri);
dguNet.plot_DGU_system(X,U, config, control_type, dguNet); % plot results

%% Add a DGU to the network: connect DGU 6 to DGU 3
filename = 'New_DGU.txt';
i = nb_subsystems+1;
[~, Vin(i), R(i),L(i),C(i), Vmax(i), Vmin(i), Imax(i), Imin(i)] = utils.importData(filename);
Rij_mat(3,6) = 1.75; Rij_mat(6,3) = Rij_mat(3,6); Rij_mat(5,6) = 0.2;
Rij_mat(6,5) = Rij_mat(5,6); % Connections from DGU 6
Vr(i) = 50.4;
Il(i) = 8.5;
dguNet2 = DGU_network(nb_subsystems+1, Rij_mat);
plot(dguNet2.graph, 'EdgeLabel', dguNet2.graph.Edges.Weight, 'Marker', 's', 'NodeColor','r', ...
      'MarkerSize', 7);
activeDGUs = 1:1:6;
for i=activeDGUs
    dguNet2 = dguNet2.initElecParam(i, Vin(i), Vr(i), Il(i), R(i), C(i), L(i),...
                              Vmax(i), Vmin(i), Imax(i), Imin(i));
end
dguNet2 = dguNet2.initDynamics(1:1:nb_subsystems+1); % recompute Dynamics with new DGU
dguNet2 = dguNet2.setConstraints(delta_config);
for i = 1:1:5
    dguNet2.Ki{i} = dguNet.Ki{i}; % copy previous passive gains
    dguNet2.Pi{i} = dguNet.Pi{i};
end
dguNet2 = PnP.redesignPhase(dguNet2, 6, "add", activeDGUs);
[x0_new, Q_Ni, Ri] = utils.tuningParam(dguNet2, delta_config);
for i=1:nb_subsystems
    x0_new{i} = X{end}(:,i); % previous states start
end

%[Xnew,Unew] = PnP.mpc_DGU_tracking(@mpc_online_2, x0_new, length_sim, dguNet2, Q_Ni, Ri);
%dguNet2.plot_DGU_system([X; Xnew],[U; Unew], config, control_type, dguNet); % plot results

%% Remove DGU 4
dguDelete = 4;
Rij_mat(dguDelete,:) = 0; Rij_mat(:,dguDelete) = 0;
activeDGUs = [1 2 3 5 6];
dguNet3 = DGU_network(6, Rij_mat);
for i=1:nb_subsystems+1
    dguNet3 = dguNet3.initElecParam(i, Vin(i), Vr(i), Il(i), R(i), C(i), L(i),...
                              Vmax(i), Vmin(i), Imax(i), Imin(i));
end
dguNet3 = dguNet3.initDynamics(1:1:nb_subsystems+1); % recompute Dynamics with new DGU
dguNet3 = dguNet3.setConstraints(delta_config);
for i = activeDGUs
    dguNet3.Ki{i} = dguNet2.Ki{i}; % copy previous passive gains
    dguNet3.Pi{i} = dguNet2.Pi{i};
end
PnP.redesignPhase(dguNet2, dguDelete, "delete", activeDGUs);
plot(dguNet3.graph);