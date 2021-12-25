%% Demonstration of a PnP operation for a DGU Network
clear
close all
% USE MOSEK as solver (ADD to path)
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'
addpath(genpath(cd));
%%
% Initialize Network: Complete configuration of the network 
utils = utilityFunctions;
PnP = SimFunctionsPnP;
config = "DISTRIBUTED"; % distributed config, for plot function
filename = 'config_DGU_1.txt';
[nb_subsystems, Vin,R,L,C, Vmax, Vmin, Imax, Imin] = utils.importData(filename);
dguNet = DGU_network(nb_subsystems); % Instantiate a DGU NETWORK class
% Set references  [V] to converge to and load current  [A]
Vr = linspace(49.975, 50.1, nb_subsystems)% references
Il = linspace(3.5, 5.5, nb_subsystems)
% set Electrical parameters and Dynamics for ALL the subsystems in the network
for i=1:nb_subsystems
    dguNet = dguNet.initElecParam(i,Vin(i), Vr(i), Il(i), R(i), C(i), L(i), ...
                                  Vmax(i), Vmin(i), Imax(i), Imin(i));
end

%% A) For now consider only 5 active DGU out of 6
activeDGU_scen1 = 1:1:5; % Initially, 5 DGUs are active
Rij_mat = zeros(nb_subsystems);
Rij_mat(1,2) = 1.75; Rij_mat(2,3) = 3.5; Rij_mat(2,4) = 1.75; 
Rij_mat(3,5) = 2.8; 
Rij_mat = Rij_mat + tril(Rij_mat',1); % Non directed graph, symmetric matrix
dguNet = dguNet.setConnectionsGraph(Rij_mat); % set links between DGU
dguNet = dguNet.setActiveDGU(activeDGU_scen1); % define which DGUs are active
figure()
plot(dguNet.NetGraph, 'EdgeLabel', dguNet.NetGraph.Edges.Weight, 'Marker', 's', 'NodeColor','r', ...
      'MarkerSize', 7);
title("Initial Scenario: 5 DGUs active");
dguNet = dguNet.initDynamics(); % initialize dynamics
% For passivity based MPC, constraints are not in  Formulation 
delta_config = false; % not in delta configuration
dguNet = dguNet.compute_Ref_Constraints(delta_config);
control_type = "MPC online";
% Use passivity to find the local passive feedback gains  and  s.t. 
dguNet = PnP.setPassiveControllers(dguNet);
% Compute Qi and Ri matrices in a distributed fashion to ensure global aymptotic
% stability
passivity = true;
[x0, Q_Ni, Ri] = utils.tuningParam(dguNet, delta_config, passivity);
% Use the tracking MPC with reconfigurable terminal ingredients  to converge to reference from the initial state
simStart = 1;
length_sim = 25;
[X, U] = PnP.mpc_DGU_tracking(@trackingMPC_reconf, x0, length_sim, dguNet, Q_Ni, Ri);
dguNet.plot_DGU_system(X,U, config, control_type, dguNet, simStart, 1:6); % plot results
% clear X U
% [X, U] = PnP.mpc_DGU_tracking(@trackingMPC_reconf_admm, x0, length_sim, dguNet, Q_Ni, Ri);
% dguNet.plot_DGU_system(X,U, config, control_type, dguNet, simStart, 1:6); % plot results

%% Test transition phase without reconfigurable terminal ingredients
delta_config = true;
use_passivity = false;
[x0, Q_Ni, Ri] = utils.tuningParam(dguNet, delta_config, use_passivity);
dguNet_delta = dguNet;
dguNet_delta = dguNet_delta.compute_Ref_Constraints(delta_config);
[dguNet_delta, Gamma_Ni, alpha_i] = offlineComputeTerminalSet(Q_Ni, Ri, dguNet_delta);
fprintf("Initial terminal set constrait alpha = %d \n", alpha_i)
alpha = alpha_i*ismember(1:6, dguNet_delta.activeDGU)';
length_sim = 30;
[Xdelt,Udelt] = PnP.mpc_sim_DGU_delta(@mpc_delta, x0, length_sim, dguNet_delta,...
                         alpha, Q_Ni, Ri, Gamma_Ni);
control_type = "MPC with offline computation of terminal ingredients";
config = "DISTRIBUTED";
simStart = 1;
dguNet.plot_DGU_system(Xdelt,Udelt, config, control_type, dguNet_delta, simStart, dguNet_delta.activeDGU)
%% B) Scenario 2: Connect DGU 6 to DGU 3
% Set all DGUs to be active. DGU 6 is now active but is not connected yet to the network
simStart2 = simStart + length_sim;
dguPos = 6;
activeDGU_scen2 = 1:1:6; % Now all the 6 DGUs are active
dguNet = dguNet.setActiveDGU(activeDGU_scen2);
dguNet2 = dguNet; % dguNet copy, with 6 active DGU but before connection
% Create connection from DGU 6 to DGU 3. For this purpose, a new instance of the network class is created with the modified structure e.g. different Laplacian matrix and  
Rij_mat(3,dguPos) = 2.75; Rij_mat(dguPos,3) = Rij_mat(3,dguPos);  % New link
dguNet2 = dguNet2.setConnectionsGraph(Rij_mat);
dguNet2 = dguNet2.initDynamics(); % recompute Dynamics (changed with integration of DGU 6)
plot(dguNet2.NetGraph, 'EdgeLabel', dguNet2.NetGraph.Edges.Weight, 'Marker', 's', 'NodeColor','r', ...
      'MarkerSize', 7);
title('Scenario 2: Connection of DGU 6 to the network')
delta_config = false;
dguNet2 = dguNet2.compute_Ref_Constraints(delta_config);

%% Redesign and Transition Phase for offline terminal ingredients
dguNet_delta = dguNet_delta.setActiveDGU(activeDGU_scen2); 
dguNet2_delta = dguNet2;
delta_config = true;
dguNet2_delta = dguNet2_delta.compute_Ref_Constraints(delta_config);
use_passivity = false;
[x0, Q_Ni, Ri, Qi] = utils.tuningParam(dguNet2_delta, delta_config, use_passivity);
[dguNet2_delta, Gamma_Ni, alpha_i] = offlineComputeTerminalSet(Q_Ni, Ri, dguNet2_delta);
fprintf("Initial terminal set constrait alpha = %d \n", alpha_i)
alpha = alpha_i*ismember(1:6, dguNet2_delta.activeDGU)';
length_sim = 30;



%% Redesign and Transition Phase for online Terminal ingredients
% Redesign Phase: Compute new  and  of neighbors set of DGU 6 (including DGU 6 itself)
dguNet2 = PnP.redesignPhase(dguNet2, dguNet2.NetGraph,dguPos, "add");
 % Re-define Q_Ni since neighbors of DGU 3 and 6 changed. Initial values for the 5 first DGUs taken from previous simulation end. 
use_passivity = true;
delta_config = false;
[x0, Q_Ni, Ri, Qi] = utils.tuningParam(dguNet2, delta_config, use_passivity); 
for i = activeDGU_scen1
    x0{i} = X{end}(:,i);   %
end
disp('x0'); celldisp(x0);
% Transition Phase: Compute steady-state value to reach to allow the plug-in of DGU 6 (PnP permitted). 
% Drive the system (the 5 initial DGU's + the 6th DGU before connection) to this steady state.
ADMM = false;
[X2_trans,U2_trans,lenSim, xs,us,alpha]= PnP.transitionPhase(x0, dguNet, dguNet2, Qi, Ri, 'reference', ADMM)
clear X2_trans U2_trans lenSim xs us alpha
ADMM = true;
[X2_trans,U2_trans,lenSim, xs,us,alpha]= PnP.transitionPhase(x0, dguNet, dguNet2, Qi, Ri, 'reference', ADMM)

% Initial states for reference MPC tracking are the states from the end of the 
% transition phase (i.e. corresponding to steady state where P&P permitted):
for i = activeDGU_scen2
    x0{i} = X2_trans{end}(:,i);   
end
lenSim2 = 25;
annot2plot.array = dguNet2.Ts*[simStart2,simStart2+lenSim];
annot2plot.text = {'Start Transition Phase', 'End Transition Phase'};
[X2, U2] = PnP.mpc_DGU_tracking(@trackingMPC_reconf_admm, x0, lenSim2, dguNet2, Q_Ni, Ri);
dguNet2.plot_DGU_system([X(1:end-1),X2_trans,X2],[U, U2_trans, U2], config, control_type, dguNet2, simStart, activeDGU_scen2, annot2plot); % plot results

%% C) 3rd Scenario: Plug out DGU 4
activeDGU_scen3 = [1 2 3 5 6]; % remove DGU 4 from active DGU list
dguNet3 = dguNet2; % copy the previous instance with 6 DGUs and create new instance for this scenario
dguNet3 = dguNet3.setActiveDGU(activeDGU_scen3);
dguDelete = 4;
Rij_mat(dguDelete,:) = 0; Rij_mat(:,dguDelete) = 0;
dguNet3 = dguNet3.setConnectionsGraph(Rij_mat);
dguNet3 = dguNet3.initDynamics();
plot(dguNet3.NetGraph, 'EdgeLabel', dguNet3.NetGraph.Edges.Weight, 'Marker', 's', 'NodeColor','r', ...
      'MarkerSize', 7);
title('Scenario where DGU 4 is to be plugged out')
dguNet3 = dguNet3.compute_Ref_Constraints(delta_config);
% Redesign Phase: Compute new  and  of neighbors set of DGU 4
dguNet3 = PnP.redesignPhase(dguNet3, dguNet2.NetGraph, dguDelete, "delete");
% Transition Phase: Take as initial state the end of simulation of scenario 2
% call again since dimension of Q_Ni change when adding/removing DGU
[~, Q_Ni, Ri, Qi] = utils.tuningParam(dguNet3, delta_config);
for i = activeDGU_scen2
    x0{i} = X2{end}(:,i);   
end
disp('x0'); celldisp(x0);
% Define new references, to see the effect of the objective function of the optimization problem:
dguNet3.Vr = linspace(49.90, 50.4, nb_subsystems);% references
dguNet3.Il = linspace(2.5, 7.5, nb_subsystems);
dguNet3 = dguNet3.compute_Ref_Constraints(delta_config);
% will keep the steady state as close as possible to the current state: for quick P&P operation
ADMM = true;
[X3_trans_,U3_trans_,lenTrans_, xs_,us_,alpha_]= PnP.transitionPhase(x0, dguNet2, dguNet3, Qi, Ri, 'current state', ADMM)
%  will keep the steady state as close as possible from the references, 
% with the goal of reducing modification to the desired system behaviour
[X3_trans,U3_trans,lenTrans, xs,us,alpha]= PnP.transitionPhase(x0, dguNet2, dguNet3, Qi, Ri, 'reference', ADMM)
% If P&P permitted, simulate normal operation of the network after plug out of DGU 4:
for i = activeDGU_scen3
    x0{i} = X3_trans{end}(:,i);   
end
lenSim3 = 20;
[X3, U3] = PnP.mpc_DGU_tracking(@trackingMPC_reconf_admm, x0, lenSim3, dguNet3, Q_Ni, Ri);
annot2plot.array = dguNet3.Ts*[11,11+lenTrans];
annot2plot.text = {'Start Transition Phase', 'End Transition Phase'};
% plot the last 11 points from the previous phase
dguNet3.plot_DGU_system([X2(end-10:end), X3_trans, X3],[U2(end-10:end), U3_trans, U3], config, ...
                control_type, dguNet3, simStart, activeDGU_scen2, annot2plot); % plot results
