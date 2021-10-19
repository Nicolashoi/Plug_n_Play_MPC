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
    %% Constants
    ma = 1; mb = ma; % mass 
    ka = 2; kb = 1.5; % spring constants
    l12 = 1; l21 = l12; % laplacian terms
    %% System dynamics
    A_1 = [0, 1; -1/ma *(ka+kb) + kb/ma*l12, 0]; B_1 = [0;1/ma];
    C_1 = [1 0]; F_1 = [0; kb/ma];
    A_2 = [0,1; -1/mb*(ka+kb), 0]; B_2 = [0; 1/mb];
    C_2 = C_1; F_2 = [0; kb/mb];
    
    % Discretization
    Ts = 10e-2;
    % first order
%     Ad_1_ = (eye(2) + Ts*A_1);
%     Bd_1_ = B_1*Ts;
%     Fd_1_ = F_1*Ts;
    %% discretization based on laplacian intercon. preserve structure
    sys1 = ss(A_1, [B_1 F_1], [], []);
    sys1_d = c2d(sys1, Ts);
    sys2 = ss(A_2, [B_2 F_2], [], []);
    sys2_d = c2d(sys2, Ts);
    
    %% Laplacian
    L = [l12, -l12; -l21, l21];
    L_tilde = L;
    C_global = [C_1, zeros(1,size(C_2,2)); zeros(1,size(C_1,2)), C_2];
%     U = L_tilde*C;
%     W = C'*L_tild;
    
    %% Non-distributed system
    A = [0, 1, 0, 0; -(ka+kb)/ma, 0, kb, 0; 0, 0, 0, 1; kb,...
         0, -(ka+kb)/mb, 0];
    B = [0, 0; 1/ma, 0; 0, 0; 0, 1/mb];
    sys = ss(A,B,C_global, []);
    sys_d = c2d(sys, Ts);
    
    %% put everything together
    param.A_1 = sys1_d.A;
    param.B_1 = sys1_d.B(:,1);
    param.F_1 = sys1_d.B(:,2);
    param.A_2 = sys2_d.A;
    param.B_2 = sys2_d.B(:,1);
    param.F_2 = sys2_d.B(:,2);
    param.C_1 = C_1;
    param.C_2 = C_2;
    param.C_global = C_global;
%     param.U = U;
%     param.W = W;
    param.L_tilde = L_tilde;
    param.global_sysd = sys_d;
    
  
end