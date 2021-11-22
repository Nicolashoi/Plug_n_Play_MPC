%% Simulate system
% Author:
%   Nicolas Hoischen
% BRIEF: % Simulate time steps of a system, given a controller in input and the
         % initial state
% Following Algorithm 2 of "Distributed Synthesis and Control of Constrained
% Linear Systems"

function [X,U, Pinf] = simulate_system(controller, x0, length_sim, simulation, param,...
                                  Q, R, Pi, Gamma_Ni, alpha)
    Pinf = [];
    switch simulation
        case "MPC"
            if param.name == "COUPLED_OSCI"
                [X,U] = mpc_sim_oscill(controller, x0, length_sim, param,...
                                          alpha, Q, R, Pi, Gamma_Ni);
            elseif param.name == "2_DGU"
                [X,U] = mpc_sim_DGU(controller, x0, length_sim, param,...
                                          alpha, Q, R, Pi, Gamma_Ni);
            elseif param.name == "DGU_delta"
                [X,U] = mpc_sim_DGU_delta(controller, x0, length_sim, param,...
                                      alpha, Q, R, Pi, Gamma_Ni);
          
            end
            
        case "PASSIVITY" % passivity only
            M = param.number_subsystem;
            Kpass = cell(1,M); D = cell(1,M); Ppass = cell(1,M); 
            Gamma_pass = cell(1,M);
            for i= 1:M
            [Kpass{i}, D{i}, Ppass{i}, Gamma_pass{i}] = controller_passivity(...
                                             param.Ai{i}, param.Bi{i},...
                                             param.Ci{i}, param.Fi{i}, ...
                                             param.L_tilde, param.global_sysd.C);
            disp("Passive controller Gain"); disp(Kpass{i});
            end
            
            if param.name == "COUPLED_OSCI"
                [X,U] = sim_oscill_distributed(x0, length_sim, Kpass,param);
            elseif param.name == "2_DGU" || param.name == "DGU_delta"
                [X,U] = sim_DGU_distributed(x0, length_sim, Kpass,param);
            end
        case "LQR"
            Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
            [Klqr,Pinf, ~] = dlqr(Ad, Bd, Q,R);
            disp("LQR controller Gain"); disp(-Klqr);
            if param.name == "COUPLED_OSCI"
                [X,U] = sim_global_coupled_osc(x0, length_sim, -Klqr,param);
            elseif param.name == "2_DGU" || param.name == "DGU_delta"
                [X,U] = sim_global_DGU(x0, length_sim,-Klqr,param);
            end
            
            
        case "MPC_2"
            M = param.number_subsystem;
            Kpass = cell(1,M); Ppass = cell(1,M); K_Ni = cell(1,M); 
            for i= 1:M
            [Kpass{i}, ~, Ppass{i}, ~] = controller_passivity(...
                                             param.Ai{i}, param.Bi{i},...
                                             param.Ci{i}, param.Fi{i}, ...
                                             param.L_tilde, param.global_sysd.C);
            disp("Passive controller Gain"); disp(Kpass{i});
            end
            for i=1:M
                neighbors = sort([i; successors(param.graph, i)]); 
                K_block = blkdiag(Kpass{neighbors}); %block matrix of all neighbors K
                K_Ni{i} = K_block(neighbors==i,:); % extract only row corresponding to subsystem i
            end
            [X,U] = mpc_DGU_tracking(controller, x0, length_sim, param,...
                                          K_Ni, Q, R, Ppass);
   
            
        otherwise
            error("Simulation not implemented");
       
    end
end




%% SIMULATION FOR OSCILLATOR
function [X,U] = sim_oscill_distributed(x0, length_sim, K,param)
    M = param.number_subsystem; % number of subsystems
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = horzcat(x0{:}); % columns are subsystem i
    for n=1:length_sim
         for i=1:M
            U{n}(:,i) = K{i}*X{n}(:,i);
            neighbors = [i, successors(param.graph, i)]; % get neighbors
            neighbors = sort(neighbors); % sorted neighbor list
            % create neighbor state vector comprising the state of subsystem i
            % and the state of it's neighbors (concatenated)
            x_Ni = reshape(X{n}(:,neighbors),[],1); 
            % Apply first input to the system and get next state
            X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
         end
    end  
end

function [X,U] = sim_global_coupled_osc(x0, length_sim, K,param)
    
    Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
    %Bd(1:2,2) = [0;0]; Bd(3:4,1) = [0;0];
    Acl = Ad+Bd*K;
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = vertcat(x0{:}); % global initial state, vector Mx1
    for n=1:length_sim
        X{n+1} = Acl*X{n};
        U{n} = K*X{n};
    end     
end


function [X, U]= mpc_sim_oscill(controller, x0, length_sim, param,...
                                          alpha, Q_Ni, Ri, Pi, Gamma_Ni)
  
        M = param.number_subsystem; % number of subsystems
        X = cell(length_sim+1,1); % state at each timestep
        U = cell(length_sim,1); % control input at each timestep
        X{1} = horzcat(x0{:}); % initial state columns are subsystem i
        N = 10; % Horizon
        clear mpc_online
        for n = 1:length_sim % loop over all subsystem
            % control input is of size nu x M
            U{n} = controller(X{n}, alpha, Q_Ni, Ri, Pi, N, param); % get first control input
            if isnan(U{n})
                error("Input to apply to controller is Nan at iteration %d",n);
            end
            for i=1:M
                neighbors = [i, successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(X{n}(:,neighbors),[],1); 
                % Apply first input to the system and get next state
                X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
                % update variable constraining terminal set of each subsystem
                alpha(i) = alpha(i) + x_Ni'*Gamma_Ni{i}*x_Ni;
            end
        end                                     
end

%% SIMULATION FOR DGU
function [X,U] = sim_DGU_distributed(x0, length_sim, K,param)
    M = param.number_subsystem; % number of subsystems
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    X{1} = horzcat(x0{:}); % columns are subsystem i
    
    if param.name == "2_DGU"
        d = cell(length_sim,1); % duty cycle of DGU's
        for n=1:length_sim
             for i=1:M
                X{n}(:,i) = X{n}(:,i) + [0;0;-param.Vr{i}];
                U{n}(:,i) = K{i}*X{n}(:,i) - param.Vr{i}/param.Vin{i} + ...
                            K{i}*[-param.Vr{i}; 0; 0];
                d{n}(:,i) = U{n}(:,i) + (param.R{i}*param.Il{i})/param.Vin{i};
                neighbors = [i, successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(X{n}(:,neighbors),[],1); 
                % Apply first input to the system and get next state
                X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
             end
        end
        U = d; % actual control input is the duty cycle
    elseif param.name == "DGU_delta"
        dX = cell(length_sim+1,1); % state at each timestep
        dU = cell(length_sim,1); % control input at each timestep
        dX{1} = X{1} - horzcat(param.Xref{:});
        
        for n=1:length_sim
             for i=1:M
                dU{n}(:,i) = K{i}*dX{n}(:,i);
                U{n} = dU{n} + param.Uref{i};
                neighbors = [i; successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(dX{n}(:,neighbors),[],1); 
                % Apply first input to the system and get next state
                dX{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*dU{n}(:,i);
                X{n+1}(:,i) = dX{n+1}(:,i) + param.Xref{i};
             end
        end  
    end 
end


function [X,U] = sim_global_DGU(x0, length_sim, K,param)
    
    Ad = param.global_sysd.A; Bd = param.global_sysd.B; 
    M = param.number_subsystem; % number of subsystems
    X = cell(length_sim+1,1); % state at each timestep
    U = cell(length_sim,1); % control input at each timestep
    d = cell(length_sim,1); % duty cycle
    X{1} = vertcat(x0{:}); % global initial state, vector Mx1
    u_to_di = cell(1,M);
    xconst = cell(1,M); uconst = cell(1,M); ufwd = cell(1,M); 
    if param.name == "2_DGU"
        for n=1:length_sim
            for i =1:M
                xconst{i} = [0;0;-param.Vr{i}];
                ufwd{i} = [-param.Vr{i}; 0; 0];
                uconst{i} = - param.Vr{i}/param.Vin{i}; 
                u_to_di{i} = (param.R{i}*param.Il{i})/param.Vin{i}; % to convert to duty cycle
            end
            X{n} = X{n} + vertcat(xconst{:});
            U{n} = K*X{n} + vertcat(uconst{:}) + K*vertcat(ufwd{:});
            X{n+1} = Ad*X{n} + Bd*U{n};
            d{n} = U{n} + vertcat(u_to_di{:});
        end     
        U = d;
     elseif param.name == "DGU_delta"
        dX = cell(length_sim+1,1); % state at each timestep
        dU = cell(length_sim,1); % control input at each timestep
        dX{1} = X{1} - vertcat(param.Xref{:});
        for n=1:length_sim
            dU{n} = K*dX{n};
            U{n} = dU{n} + vertcat(param.Uref{:});
            dX{n+1}= Ad*dX{n}+ Bd*dU{n};
            X{n+1} = dX{n+1} + vertcat(param.Xref{:});
        end
    end
end

function [X,U] = mpc_sim_DGU_delta(controller, x0, length_sim, param,...
                                          alpha, Q_Ni, Ri, Pi, Gamma_Ni)
                                      
        M = param.number_subsystem; % number of subsystems
        dX = cell(length_sim+1,1); % state at each timestep
        dU = cell(length_sim,1); % control input at each timestep
        X = cell(length_sim+1,1); % state at each timestep
        U = cell(length_sim,1); % control input at each timestep
        X{1} = horzcat(x0{:}); % initial state columns are subsystem i    
        dX{1} = X{1} - horzcat(param.Xref{:});
        N = 10; % Horizon
        clear mpc_online
        for n = 1:length_sim % loop over all subsystem
            % control input is of size nu x M
            dU{n} = controller(dX{n}, alpha, Q_Ni, Ri, Pi, N, param); % get first control input
            U{n} = dU{n} + horzcat(param.Uref{:});
            if isnan(dU{n})
                error("Input to apply to controller is Nan at iteration %d",n);
            end
            for i=1:M
                neighbors = [i; successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(dX{n}(:,neighbors),[],1); 
                % Apply first input to the system and get next state
                dX{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*dU{n}(:,i);
                X{n+1}(:,i) = dX{n+1}(:,i) + param.Xref{i};
                % update variable constraining terminal set of each subsystem
                alpha(i) = alpha(i) + x_Ni'*Gamma_Ni{i}*x_Ni;
            end
        end                                                                             
end

function [X,U] = mpc_DGU_tracking(controller, x0, length_sim, param,...
                                          K, Q_Ni, Ri, Pi)
                                      
        M = param.number_subsystem; % number of subsystems
        X = cell(length_sim+1,1); % state at each timestep
        U = cell(length_sim,1); % control input at each timestep
        X{1} = horzcat(x0{:}); % initial state columns are subsystem i    
   
        N = 10; % Horizon
        clear mpc_online_2
        for n = 1:length_sim % loop over all subsystem
            % control input is of size nu x M 
            U{n} = controller(X{n}, K, Q_Ni, Ri, Pi, N, param); % get first control input
            if isnan(U{n})
                error("Input to apply to controller is Nan at iteration %d",n);
            end
            for i=1:M
                neighbors = [i; successors(param.graph, i)]; % get neighbors
                neighbors = sort(neighbors); % sorted neighbor list
                % create neighbor state vector comprising the state of subsystem i
                % and the state of it's neighbors (concatenated)
                x_Ni = reshape(X{n}(:,neighbors),[],1); 
                % Apply first input to the system and get next state
                X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i);
                %X{n}(2,i) = X{n}(2,i) + param.Il(i);
            end
        end  
       % X{end}(2,:) = X{n}(2,:) + param.Il(:);
end

% function [X, U]= mpc_sim_DGU(controller, x0, length_sim, param,...
%                                           alpha, Q_Ni, Ri, Pi, Gamma_Ni)
%   
%         M = param.number_subsystem; % number of subsystems
%         X = cell(length_sim+1,1); % state at each timestep
%         U = cell(length_sim,1); % control input at each timestep
%         X{1} = horzcat(x0{:}); % initial state columns are subsystem i
%         sfi = cell(1,M);
%         Kpass = cell(1,M);
%         for i=1:M
%             Kpass{i} = controller_passivity(param.Ai{i}, param.Bi{i},...
%                                              param.Ci{i}, param.Fi{i}, ...
%                                              param.L_tilde, param.global_sysd.C);
%             sfi{i} = [0;0;-param.Vr{i}];
%         end
%         sf = horzcat(sfi{:});
%         X{1} = X{1} + sf;
%         N = 10; % Horizon
%         clear mpc_online
%         for n = 1:length_sim % loop over all subsystem
%             control input is of size nu x M
%             U{n} = controller(X{n}, alpha, Q_Ni, Ri, Pi, N, param); % get first control input
%             if isnan(U{n})
%                 error("Input to apply to controller is Nan at iteration %d",n);
%             end
%             for i=1:M
%                 neighbors = [i, successors(param.graph, i)]; % get neighbors
%                 neighbors = sort(neighbors); % sorted neighbor list
%                 create neighbor state vector comprising the state of subsystem i
%                 and the state of it's neighbors (concatenated)
%                 x_Ni = reshape(X{n}(:,neighbors),[],1); 
%                 Apply first input to the system and get next state
%                 U{n}(:,i) = - param.Vr{i}/param.Vin{i} + Kpass{i}*[-param.Vr{i};0;0];
%                 X{n+1}(:, i) = param.A_Ni{i}*x_Ni + param.Bi{i}*U{n}(:,i) - sfi{i};
%                 update variable constraining terminal set of each subsystem
%                 alpha(i) = alpha(i) + x_Ni'*Gamma_Ni{i}*x_Ni;
%             end
%         end                                     
% end
