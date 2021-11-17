%% Parameter file for a Network of 2 DGUs %%
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   2 DGU connected by a resistor

function param = param_DGU_delta
    %% Converter Parameters
    DGU_PARAM_ELEC
    %% Subystem dynamics
    Ts = 1e-5; % sampling time
    Ac{1} = [0, 1/Ct{1}; -1/Li{1}, -R{1}/Li{1}];
    Bc{1} = [0; Vin{1}/Li{1}];
    Fc{1} = [1/Ct{1}; 0];
    Cc{1} = [1, 0];
    Ac{2} = [0, 1/Ct{2}; -1/Li{2}, -R{2}/Li{2}];
    Bc{2} = [0; Vin{2}/Li{2}];
    Fc{2} = [1/Ct{2}; 0];
    Cc{2} = Cc{1};
    
    utils = utilityFunctions;
   [di_ref, Iti_ref] = utils.compute_ref();
    %% Global system dynamics
    A = blkdiag(Ac{1}, Ac{2});
    % add interconnections
    a12 = Agraph(1,2); a21 = Agraph(2,1);
    A = A + [-a12*Fc{1}*Cc{1}, a12*Fc{1}*Cc{2}; a21*Fc{2}*Cc{1},-a21*Fc{2}*Cc{2}]; 
    B = blkdiag(Bc{1},Bc{2});
    C = blkdiag(Cc{1},Cc{2});
    sys = ss(A,B,C,[]);
    sys_d = c2d(sys, Ts); % Exact discretization of the global system
    %% Laplacian
    
    L_tilde = kron(L, eye(size(Cc{1},1)));
    
    %% Discretization of subsystems
    for i = nb_subsystems:-1:1 % in descending order so that Matlab allocates free space
        subsys = ss(Ac{i}, [Bc{i} Fc{i}], [], []);
        subsys_d{i} = c2d(subsys, Ts);
        Ai{i} = subsys_d{i}.A;
        Bi{i} = subsys_d{i}.B(:,1);
        Fi{i} = subsys_d{i}.B(:,2);
        Ci{i} = Cc{i};
    end
    
    A_Ni = utils.change_system_representation(Ai,Fi,Ci,Agraph);
    %% Parameters for offline algorithm 1
    U{1} = [1 0 0 0; 0 1 0 0]; 
    U{2} = [0 0 1 0; 0 0 0 1];
    W{1} = eye(4); W{2} = eye(4);
    
    %% constraints now in delta formulation
    Xref = cell(1,nb_subsystems);
    Uref = cell(1,nb_subsystems);
    Gx_i = cell(1,nb_subsystems);
    Gu_i = cell(1,nb_subsystems);
    fx_i = cell(1,nb_subsystems);
    fu_i = cell(1,nb_subsystems);
    for i= 1:nb_subsystems
        Xref{i} = [Vr{i}; Iti_ref{i}];
        Uref{i} = di_ref{i};
        Gx_i{i}= [1,0;-1 0; 0 1; 0 -1];
        Gu_i{i} = [1;-1];
        fx_i{i} = [60; -40; 10; 0] + [-Xref{i}; Xref{i}];
        fu_i{i} = [1;0] + [-Uref{i}; Uref{i}];
    end 
   
    Gx = blkdiag(Gx_i{:});
    fx = [fx_i{:}];
    Gu = blkdiag(Gu_i{:});
    fu = [fu_i{:}];
    
    graph = digraph([1 2], [2 1]); % or digraph(Agraph)
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
    param.size_subsystem = 2;
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