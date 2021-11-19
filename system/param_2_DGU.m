%% Parameter file for a Network of 2 DGUs %%
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   2 DGU connected by a resistor

function param = param_2_DGU
    %% Converter Parameters
    R{1} = 3e-3; L1 = 92e-6; Cap1 = 75e-6;
    R{2} = 1.5e-3; L2 = 102e-6; Cap2 = 50e-6;
    R12 = 1.75; 
    Vin{1} = 100; Vin{2} = 100;
    Vr{1} = 50; Il{1} = 5; 
    Vr{2} = 50; Il{2} = 5;
    nb_subsystems = 2; 
    a11 = 0; a12 = 1/R12; a21 = a12; a22 = 0;
    %% Subystem dynamics
    Ts = 1e-5; % sampling time
    Ac{1} = [0, 1/Cap1, 0; -1/L1, -R{1}/L1, 0; 1/Ts, 0, 0];
    Bc{1} = [0; Vin{1}/L1; 0];
    Fc{1} = [1/Cap1; 0; 0];
    Cc{1} = [1, 0, 0];
    Ac{2} = [0, 1/Cap2, 0; -1/L2, -R{2}/L2, 0; 1/Ts, 0, 0];
    Bc{2} = [0; Vin{2}/L2; 0];
    Fc{2} = [1/Cap2; 0; 0];
    Cc{2} = Cc{1};
    
    %% Global system dynamics
    A = blkdiag(Ac{1}, Ac{2});
    % add interconnections
    A = A + [-a12*Fc{1}*Cc{1}, a12*Fc{1}*Cc{2}; a21*Fc{2}*Cc{1},-a21*Fc{2}*Cc{2}]; 
    B = blkdiag(Bc{1},Bc{2});
    C = blkdiag(Cc{1},Cc{2});
    sys = ss(A,B,C,[]);
    sys_d = c2d(sys, Ts); % Exact discretization of the global system
    %% Laplacian
    L = [a11+a12,-a12; -a21, a21 + a22];
    L_tilde = kron(L, eye(size(Cc{1},1)));
    Agraph = [0, a12; a21, 0]; % graph matrix
    %% Discretization of subsystems
    for i = nb_subsystems:-1:1 % in descending order so that Matlab allocates free space
        subsys = ss(Ac{i}, [Bc{i} Fc{i}], [], []);
        subsys_d{i} = c2d(subsys, Ts);
        Ai{i} = subsys_d{i}.A;
        Bi{i} = subsys_d{i}.B(:,1);
        Fi{i} = subsys_d{i}.B(:,2);
        Ci{i} = Cc{i};
    end
    utils = utilityFunctions;
    A_Ni = utils.change_system_representation(Ai,Fi,Ci,Agraph);
    %% Parameters for offline algorithm 1
    U{1} = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0]; 
    U{2} = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
    W{1} = eye(6); W{2} = eye(6);
    
    %% constraints
    Gx_i{1}= [eye(3); -eye(3)]; Gx_i{2} = Gx_i{1};
    fx_i{1} = [80;10;1000000; -30; 0; 100000]; fx_i{2} = fx_i{1};
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