%% Parameter file for a coupled oscillator %%
% Author:
%   Nicolas Hoischen
% BRIEF: 
%   Parameter definition for a simple coupled oscillator model
%   Two masses ma and mb
%   each connected to the wall with a spring of constant ka
%   interconnection of masses through spring of constant kb
%   xa = displacement of first mass from equilibrium
%   xb = displacement of second mass from equilibrium

function param = param_coupled_oscillator
    %% Constant parameters
    m1 = 1; m2 = m1; % mass 
    ka = 2; kb = 1.5; % spring constants
    a12 = 1; a21 = a12; % laplacian terms
    nb_subsystems = 2;   
    Ts = 10e-2; % Discretization
    %% System dynamics
    Ac{1} = [0, 1; -ka/m1, 0]; Bc{1} = [0;1/m1];
    Cc{1} = [1 0]; Fc{1} = [0; kb/m1];
    Ac{2} = [0,1; -ka/m2, 0]; Bc{2} = [0; 1/m2];
    Cc{2} = Cc{1}; Fc{2} = [0; kb/m2];
    
    %% discretization based on laplacian intercon. to preserve structure
    for i = nb_subsystems:-1:1 % in descending order so that Matlab allocates free space
        sys = ss(Ac{i}, [Bc{i} Fc{i}], [], []);
        sys_d{i} = c2d(sys, Ts);
        Ai{i} = sys_d{i}.A;
        Bi{i} = sys_d{i}.B(:,1);
        Fi{i} = sys_d{i}.B(:,2);
        Ci{i} = Cc{i};
    end
     %% Non-distributed system
    A = [0, 1, 0, 0; -(ka+kb)/m1, 0, kb, 0; 0, 0, 0, 1; kb,...
         0, -(ka+kb)/m2, 0];
    B = [0, 0; 1/m1, 0; 0, 0; 0, 1/m2];
    C = blkdiag(Cc{1}, Cc{2});
    sys = ss(A,B,C,[]);
    sys_d = c2d(sys, Ts); % Exact discretization of the global system
    A_Ni{1} = sys_d.A(1:2,:);
    A_Ni{2} = sys_d.A(3:4,:);
    %% Laplacian
    Agraph = [0, a12; a21, 0]; % graph matrix
    L = [a12, -a12; -a21, a21];
    L_tilde = L;
    
    %% Parameters for offline algorithm 1
    U{1} = [1 0 0 0; 0 1 0 0]; U{2} = [0 0 1 0; 0 0 0 1];
    W{1} = eye(4); W{2} = eye(4);
    
    % constraints
    Gx_i{1}= [1,0; -1,0]; Gx_i{2} = Gx_i{1};
    fx_i{1} = [1;1]; fx_i{2} = fx_i{1};
    Gu_i{1} = [1;-1]; Gu_i{2} = Gu_i{1};
    fu_i{1} = [1;1]; fu_i{2} = fu_i{1};
    Gx = blkdiag(Gx_i{1}, Gx_i{2});
    fx = [fx_i{1}; fx_i{2}];
    Gu = blkdiag(Gu_i{1}, Gu_i{2});
    fu = [fu_i{1}; fu_i{2}];
    
    graph = digraph([1 2], [2 1]);
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
    param.name = "COUPLED_OSCI";
    param.A_Ni = A_Ni;
    param.graph = graph;
    
  
end