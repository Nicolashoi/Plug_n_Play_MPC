%% Parameter file for a Network of 2 DGUs %%
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   2 DGU connected by a resistor

function param = param_DGU_delta
    %% Converter Parameters
    load('system/DGU_electrical_param.mat')
    % references
    Vr = linspace(49.95, 50.05, 6);
    Il = linspace(2.5, 7.5, 6);
    %% Subystem dynamics
    Ts = 1e-5; % sampling time
    for i = 1:nb_subsystems
        Ac{i} = [0, 1/Ct(i); -1/Li(i), -R(i)/Li(i)];
        Bc{i} = [0; Vin(i)/Li(i)];
        Fc{i} = [1/Ct(i); 0];
        Cc{i} = [1, 0];
    end
    % Discretization of subsystems
    for i = nb_subsystems:-1:1 % in descending order so that Matlab allocates free space
        subsys = ss(Ac{i}, [Bc{i} Fc{i}], [], []);
        subsys_d{i} = c2d(subsys, Ts);
        Ai{i} = subsys_d{i}.A;
        Bi{i} = subsys_d{i}.B(:,1);
        Fi{i} = subsys_d{i}.B(:,2);
        Ci{i} = Cc{i};
    end
    utils = utilityFunctions;
    % A_Ni matrices in discrete time 
    A_Ni = utils.change_system_representation(Ai,Fi,Ci,Agraph);
    % Compute references
   [di_ref, Iti_ref] = utils.compute_ref(nb_subsystems,Agraph, Vr, Il, Rij, R, Vin);
    %% Laplacian
    %L_tilde1 = kron(Agraph, eye(size(Cc{1},1)));
    L_tilde = L; % laplacian
    %% Global system dynamics
    A = blkdiag(Ac{:});
    B = blkdiag(Bc{:});
    C = blkdiag(Cc{:});
    F = blkdiag(Fi{:});
    A = A + F*L_tilde*C;
%     A_Ni_c = utils.change_system_representation(Ac,Fc,Cc,Agraph);
%     %paddzeros = max(cellfun(@length, A_Ni_c));
%     size_xglobal = nb_subsystems * size_subsystem;
%     % padd zeros to concatenate A_Ni_c vertically
%     A_Ni_c = cellfun(@(x)[x,zeros(2, size_xglobal-size(x,2))],A_Ni_c, ...
%                          'UniformOutput',false);
%     A = vertcat(A_Ni_c{:});
%     
    sys = ss(A,B,C,[]);
    sys_d = c2d(sys, Ts); % Exact discretization of the global system
  
    %% Parameters for offline algorithm 1
    %% CHANGE THIS TO 6 SUBSYSTEMS -> MOVE TO MAT FILE ?
    U = cell(1,nb_subsystems);
    W = cell(1,nb_subsystems);
    mi = size_subsystem;
    N = nb_subsystems* mi;
    state_i = eye(2);
    graph = digraph(Agraph); % or digraph(Agraph)
    for i = 1:nb_subsystems
        U{i} = [zeros(mi,mi*(i-1)), state_i, zeros(mi,N-mi*i)];
        out_neighbors = sort([i;successors(graph, i)]); % neighbor states of i
        W{i} = zeros(nb_subsystems);
        % put a 1 in diagonal of neighbor state
        for j=1:length(out_neighbors)
            W{i}(out_neighbors(j),out_neighbors(j)) = 1;
        end
        W{i}= kron(W{i}, state_i); % each subsystem has size mi
        W{i}(all(W{i}==0,2),:)=[]; % remove all zero rows
    end
   
    %% constraints now in delta formulation
    Xref = cell(1,nb_subsystems);
    Uref = cell(1,nb_subsystems);
    Gx_i = cell(1,nb_subsystems);
    Gu_i = cell(1,nb_subsystems);
    fx_i = cell(1,nb_subsystems);
    fu_i = cell(1,nb_subsystems);
    for i= 1:nb_subsystems
        Xref{i} = [Vr(i); Iti_ref(i)];
        Uref{i} = di_ref(i);
        Gx_i{i}= [eye(2); -eye(2)];
        Gu_i{i} = [1;-1];
        fx_i{i} = [52; 10; -49; 0] + [-Xref{i}; Xref{i}];
        fu_i{i} = [1;0] + [-Uref{i}; Uref{i}];
    end 
   
    Gx = blkdiag(Gx_i{:});
    fx = [fx_i{:}];
    Gu = blkdiag(Gu_i{:});
    fu = [fu_i{:}];
    
    
    %% put everything together
    param.Ts= Ts;
    param.Ai = Ai;
    param.Bi = Bi;
    param.Fi = Fi;
    param.Ci = Ci;
    param.L_tilde = L_tilde;
    param.global_sysd = sys_d;
    param.U = U;
    param.W = W;
    param.Gx_i = Gx_i;
    param.Gu_i = Gu_i;
    param.fx_i = fx_i;
    param.fu_i = fu_i;
    param.Gx = Gx;
    param.Gu = Gu;
    param.fx = fx;
    param.fu = fu;
    param.Agraph = Agraph;
    param.size_subsystem = size_subsystem;
    param.number_subsystem = nb_subsystems;
    param.Vr = Vr;
    param.Il = Il;
    param.R = R;
    param.Vin= Vin;
    param.name = "DGU_delta";
    param.A_Ni = A_Ni;
    param.graph = graph;
    param.Xref = Xref;
    param.Uref = Uref;
end