%%
clear
% USE MOSEK as solver (ADD to path)
addpath 'C:\Program Files\Mosek\9.3\toolbox\R2015aom'
addpath(genpath(cd));

%% PART 1: PASSIVITY VS LQR
length_sim = 50;
simStart = 1;
utils = utilityFunctions; % utility function class
sim = SimFunctionsPnP; % simulation function class
filename = 'config_DGU_1.txt'; % Initial DGU Network
[nb_subsystems, Vin,R,L,C, Vmax, Vmin, Imax, Imin] = utils.importData(filename);
activeDGU = 1:nb_subsystems; % all DGU are active
Rij_mat = zeros(nb_subsystems);
Rij_mat(1,2) = 1.75; Rij_mat(2,3) = 3.5; Rij_mat(2,4) = 1.75; 
Rij_mat(3,5) = 2.8;
Rij_mat = Rij_mat + tril(Rij_mat',1); % symmetric matrix for undirected graph
dguNet = DGU_network(nb_subsystems);
dguNet = dguNet.setConnectionsGraph(Rij_mat);
dguNet = dguNet.setActiveDGU(activeDGU);
Vr = linspace(49.95, 50.2, nb_subsystems); % references
Il = linspace(2.5, 7.5, nb_subsystems);
% set Electrical parameters
for i=1:nb_subsystems
    dguNet = dguNet.initElecParam(i,Vin(i), Vr(i), Il(i), R(i), C(i), L(i), ...
                                  Vmax(i), Vmin(i), Imax(i), Imin(i));
end
dguNet = dguNet.initDynamics();
R36 = 0.005:0.05:5; % no QI, Ri found for R36 < 0.6 for passivity to guarante asympt. stability
Qi = cell(1,length(R36)); Ri = cell(1,length(R36)); lambda = zeros(1, length(R36));
costPassivity = zeros(1,length(R36)); costLQR = zeros(1,length(R36));
feasiblePass = zeros(length(R36),1);
feasibleLyap = zeros(length(R36),1);
for i = 1:length(R36)
    Rij_mat(3,6) = R36(i); Rij_mat(6,3) = R36(i);
    dguNet = dguNet.setConnectionsGraph(Rij_mat);
    dguNet = dguNet.initDynamics();
    delta_config = true; % Delta-Formulation of DGU Network
    dguNet = dguNet.compute_Ref_Constraints(delta_config); % constraints in delta formulation
    control_type = "PASSIVITY";
    config = "DISTRIBUTED";
    [~, Q_Ni, Ri, Qi] = utils.tuningParam(dguNet, delta_config, false);
    [~, ~, ~,feasibleLyap(i)] = offlineComputeTerminalSet(Q_Ni, Ri, dguNet);
    [dguNet, lambda(i), feasiblePass(i)] = sim.setPassiveControllers(dguNet);
    
    fprintf("%d feasible states for passivity out of %d tested connection strengths \n ", sum(feasiblePass), length(R36));
    fprintf("%d feasible states for lyapunov out of %d tested connection strengths \n", sum(feasibleLyap), length(R36));
%     passivity = true; % set passivity boolean to true
%     [x0, ~, Ri, Qi] = utils.tuningParam(dguNet, delta_config, passivity); % initial parameters definition
%     [X,U] = sim.sim_DGU_distributed(x0, length_sim, dguNet, dguNet.Ki);
%     error_pass(i) = utils.tracking_error(X,U,config,dguNet,activeDGU);
%     % LQR
%     config = "GENERAL";
%     control_type = 'LQR';
%     Ad = dguNet.global_sysd.A; Bd = dguNet.global_sysd.B; 
%     Q = eye(size(Ad,1)); R = eye(size(Bd,2));
%     disp("LQR controller Gain");
%     [Klqr, Pinf, ~] = dlqr(Ad, Bd, Q,R);
%     [Xlqr,Ulqr] = sim.sim_global_DGU(x0, length_sim, dguNet, -Klqr);
%     error_lqr(i) = utils.tracking_error(Xlqr,Ulqr,config,dguNet,activeDGU);
%     % Cost
%     Qpass = blkdiag(Qi{:}); Rpass = blkdiag(Ri{:});
%     costPassivity(i) = utils.compute_QR_cost(X,U,Q,R,...
%                                         "DISTRIBUTED", dguNet.Xref , dguNet.Uref);
%     costLQR(i) = utils.compute_QR_cost(Xlqr,Ulqr,Q,R,"GENERAL",dguNet.Xref,dguNet.Uref);
end
% suboptIdx = (error_pass-error_lqr)./error_lqr;
% meanSubopt = mean(suboptIdx);
% stdSubopt = std(suboptIdx);
% meanCostPassivity = mean(costPassivity); stdCostPassivity = std(costPassivity);
% meanCostLQR = mean(costLQR); stdCostLQR = std(costLQR);
disp("Min eigenvalue dissipation rate"); disp(lambda);
%% PART 2: Transition Phase
clear
utils = utilityFunctions;
PnP = SimFunctionsPnP;
config = "DISTRIBUTED"; % distributed config, for plot function
filename = 'config_DGU_1.txt';
[nb_subsystems, Vin,R,L,C, Vmax, Vmin, Imax, Imin] = utils.importData(filename);
dguNet = DGU_network(nb_subsystems); % Instantiate a DGU NETWORK class
Vr = linspace(49.95, 50.05, nb_subsystems);% references

costFunction = 'reference';
MC_iter = 10;
CostRelativeDiff = zeros(MC_iter,1);
[xs, us, xs_delt, us_delt] = deal(cell(1,MC_iter), cell(1,MC_iter), cell(1,MC_iter),...
                                  cell(1,MC_iter));
[SolverTime,SolverTimeDelta] = deal(zeros(MC_iter,1), zeros(MC_iter,1));
maxVr = 0.5; minVr = -0.5; 
maxIl = 7.5; minIl = 2.5;

Il = linspace(2.5, 7.5, nb_subsystems);

for k=1:MC_iter
    Vdrop = (maxVr - minVr)*rand(1,dguNet.nb_subsystems)+ minVr;
    %Il = (maxIl - minIl)*rand(1,dguNet.nb_subsystems) + minIl;
    x0 = cell(1, dguNet.nb_subsystems); x0_delta = cell(1,dguNet.nb_subsystems);

    % set Electrical parameters and Dynamics for ALL the subsystems in the network
    for i=1:nb_subsystems
        dguNet = dguNet.initElecParam(i,Vin(i), Vr(i), Il(i), R(i), C(i), L(i), ...
                                      Vmax(i), Vmin(i), Imax(i), Imin(i));
    end
    for i = 1:dguNet.nb_subsystems
        x0_delta{i} = [50+Vdrop(i);5]; % initial condition in normal coordinates
        x0{i} = [50+Vdrop(i);5-dguNet.Il(i)]; % second state is Ii - Il
    end
    %--------------------- 5 DGU ACTIVE ------------------------------------------%
    activeDGU_scen1 = 1:1:5; % Initially, 5 DGUs are active
    Rij_mat = zeros(nb_subsystems);
    Rij_mat(1,2) = 1.75; Rij_mat(2,3) = 3.5; Rij_mat(2,4) = 1.75; 
    Rij_mat(3,5) = 2.8; 
    Rij_mat = Rij_mat + tril(Rij_mat',1); % Non directed graph, symmetric matrix
    dguNet = dguNet.setConnectionsGraph(Rij_mat); % set links between DGU
    dguNet = dguNet.setActiveDGU(activeDGU_scen1); % define which DGUs are active
    dguNet = dguNet.initDynamics(); % initialize dynamics
    % set Passivity Gains
    delta_config = false; % not in delta configuration
    dguNet = dguNet.compute_Ref_Constraints(delta_config);
    dguNet = PnP.setPassiveControllers(dguNet);
    passivity = true;
    [~, ~, Ri, Qi] = utils.tuningParam(dguNet, delta_config, passivity);
    % Delta with offline computation
    delta_config = true;
    dguNet_delta = dguNet;
    dguNet_delta = dguNet_delta.compute_Ref_Constraints(delta_config);
    %------- Connect DGU 6 to DGU 3 -----------------------------------------------%
    dguPos = 6;
    activeDGU_scen2 = 1:1:6; % Now all the 6 DGUs are active
    dguNet = dguNet.setActiveDGU(activeDGU_scen2);
    dguNet2 = dguNet; % dguNet copy, with 6 active DGU but before connection
    Rij_mat(3,dguPos) = 2.75; Rij_mat(dguPos,3) = Rij_mat(3,dguPos);  % New link
    dguNet2 = dguNet2.setConnectionsGraph(Rij_mat);

    dguNet2 = dguNet2.initDynamics(); % recompute Dynamics (changed with integration of DGU 6)
    delta_config = false;
    dguNet2 = dguNet2.compute_Ref_Constraints(delta_config);
    [dguNet2, Qi, Ri, Q_Ni] = PnP.redesignPhase(dguNet2, dguNet2.NetGraph,dguPos, "add", Qi, Ri);

    dguNet_delta = dguNet_delta.setActiveDGU(activeDGU_scen2); 
    dguNet2_delta = dguNet2;
    delta_config = true;
    dguNet2_delta = dguNet2_delta.compute_Ref_Constraints(delta_config);
    use_passivity = false;
    [~, Q_Ni_delt, Ri_delt, Qi_delt] = utils.tuningParam(dguNet2_delta, delta_config, use_passivity);
    [dguNet2_delta, Gamma_Ni, alpha_i] = offlineComputeTerminalSet(Q_Ni_delt, Ri_delt, dguNet2_delta);
    fprintf("Initial terminal set constrait alpha = %d \n", alpha_i)
    alpha_delt = alpha_i*ismember(1:6, dguNet2_delta.activeDGU)';

    ADMM = true; regulation = false;
    [~, ~, ~, xs{k},us{k},alpha, SolverTime(k)] = PnP.transitionPhase(x0, dguNet, dguNet2, Qi, Ri, costFunction, ADMM, regulation);
    [~, ~, ~, xs_delt{k},us_delt{k}, SolverTimeDelta(k)] = PnP.transitionPhaseDeltaADMM(x0_delta, dguNet_delta,...
                                      dguNet2_delta, Qi_delt, Ri_delt,costFunction, alpha_delt, regulation);
    %xs = xs(2,:) + dguNet2.Il; % add load currents
    if strcmp(costFunction, 'reference')
        costOnlineTransition =  norm(xs{k}-horzcat(dguNet2.Xref{:}),2);
        costDeltaTransition =  norm(xs_delt{k}-horzcat(dguNet2_delta.Xref{:}),2);
    elseif strcmp(costFunction, 'current state')
        costOnlineTransition =  norm(xs{k}-horzcat(x0{:}),2);
        costDeltaTransition =  norm(xs_delt{k}-horzcat(x0_delta{:}),2);
    else
        disp("cost function not well defined");
    end
    CostRelativeDiff(k) = (costOnlineTransition - costDeltaTransition)./min(costOnlineTransition, costDeltaTransition)*100;
end
fprintf("Mean solver time for online terminal ingredients = %d with std = %d \n", ...
        mean(SolverTime), std(SolverTime));
fprintf("Mean solver time for offline terminal ingredients = %d with std = %d", ...
        mean(SolverTimeDelta), std(SolverTimeDelta));