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
    U{1} = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0]; 
    U{2} = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
    W{1} = eye(6); W{2} = eye(6);
    
    %% constraints
    Gx_i{1}= [1,0, 0; -1,0, 0; 0, 1, 0; 0, -1, 0]; Gx_i{2} = Gx_i{1};
    fx_i{1} = [500;500; 200; 200]; fx_i{2} = fx_i{1};
    Gu_i{1} = [1;-1]; Gu_i{2} = Gu_i{1};
    fu_i{1} = [1000;1000]; fu_i{2} = fu_i{1};
    Gx = blkdiag(Gx_i{1}, Gx_i{2});
    fx = [fx_i{1}; fx_i{2}];
    Gu = blkdiag(Gu_i{1}, Gu_i{2});
    fu = [fu_i{1}; fu_i{2}];
    
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
    param.size_subsystem = 3;
    param.number_subsystem = nb_subsystems;
    param.Vr = Vr;
    param.Il = Il;
    param.R = R;
    param.Vin= Vin;
    param.name = "2_DGU";
    param.A_Ni = A_Ni;
    param.graph = graph;
    
end