%% Generic class to define a DGU_network
% Author:
%   Nicolas Hoischen
% BRIEF:
    % Class to define a DGU network with different methods to create new
    % references, new connections etc..
    % Automatically computes the dynamics and references given electrical
    % parameters as input
    % Also contains a method to plot the DGU state & duty cycle evolution 
    % Handles both normal formulation and delta formulation of the dynamics

classdef DGU_network
    properties
        
       %% Electrical Parameters
       Ts = 1e-5;
       Vin 
       Vr
       Il 
       Ri 
       Ct 
       Li 
       %% Dynamics
       Ai
       Bi
       Fi
       Ci
       A_Ni
       global_sysd % useful to compare to LQR 
       
       %% Controller Parameters
       Ki
       K_Ni
       Pi
       %% Others
       Agraph
       NetGraph
       Rij_mat
       L_tilde
       Iti_ref
       di_ref
       nb_subsystems
       ni = 2 % size of subsystem state
       nu = 1
       delta_config
       %% Constraints
       Vmax
       Vmin
       Imax
       Imin
       Gx_i
       Gu_i
       fx_i
       fu_i
       Gx_Ni
       fx_Ni
       Xref
       Uref
      
       %% param for state selection
       U
       W
       Wij
       activeDGU
    end
    properties (Access = private)
       Ac
       Bc
       Fc
       Cc 
    end
    
    methods
        % constructor
        function obj = DGU_network(nb_subsystems)
            if nargin == 1
                obj.nb_subsystems = nb_subsystems;
            end
        end
        % define which DGUs are active
        function obj = setActiveDGU(obj, activeDGU)
            obj.activeDGU = activeDGU;
        end
        % Set the interconnection strengths in the network, Laplacian & Graph
        % Matrix
        function obj = setConnectionsGraph(obj, Rij_mat)
            obj.Rij_mat = Rij_mat;
            obj.Agraph = Rij_mat;
            nonzeroIdx = Rij_mat ~= 0;
            obj.Agraph(nonzeroIdx) = 1./obj.Agraph(nonzeroIdx);  
            obj.NetGraph = graph(obj.Agraph);
            obj = obj.setSelectionMatrices(obj);
            obj.L_tilde = diag(sum(obj.Agraph))-obj.Agraph;
        end
        
        % initialize electrical parameters
        function obj = initElecParam(obj, idx_DGU, Vin, Vr, Il, Ri, Ct, Li, ...
                                     Vmax, Vmin, Imax, Imin)
            obj.Vin(idx_DGU) = Vin;
            obj.Vr(idx_DGU) = Vr;
            obj.Il(idx_DGU) = Il;
            obj.Ri(idx_DGU) = Ri;
            obj.Ct(idx_DGU) = Ct;
            obj.Li(idx_DGU) = Li;
            obj.Vmax(idx_DGU) = Vmax;
            obj.Vmin(idx_DGU) = Vmin;
            obj.Imax(idx_DGU) = Imax;
            obj.Imin(idx_DGU) = Imin;
        end
        
        % initialize DGU dynamics
        function obj = initDynamics(obj)
            for i = 1:obj.nb_subsystems  
                obj.Ac{i} = [0, 1/obj.Ct(i); -1/obj.Li(i), -obj.Ri(i)/obj.Li(i)];
                obj.Bc{i} = [0; obj.Vin(i)/obj.Li(i)];
                obj.Fc{i} = [1/obj.Ct(i); 0];
                obj.Cc{i} = [1, 0];  
                subsys = ss(obj.Ac{i}, [obj.Bc{i} obj.Fc{i}], [], []);
                subsys_d = c2d(subsys, obj.Ts);
                obj.Ai{i} = subsys_d.A;
                obj.Bi{i} = subsys_d.B(:,1);
                obj.Fi{i} = subsys_d.B(:,2);
                obj.Ci{i} = obj.Cc{i};     
            end
                A = blkdiag(obj.Ac{:});
                B = blkdiag(obj.Bc{:});
                C = blkdiag(obj.Cc{:});% define the global output matrix needed for passivity
                F = blkdiag(obj.Fc{:});
                A = A - F*obj.L_tilde*C;
                sys = ss(A,B,C,[]);
                obj.global_sysd = c2d(sys, obj.Ts); % Exact discretization of the global system
                
            % function to change representation to A_Ni
            obj.A_Ni = change_system_representation(obj.Ai, obj.Fi, obj.Ci, obj.Agraph);
            % function compute references here
            [obj.di_ref, obj.Iti_ref] = compute_ref(obj.nb_subsystems, obj.Agraph,...
                                                    obj.Vr, obj.Il, obj.Rij_mat, obj.Ri, obj.Vin); 
        end
        % computes the references depending if the DGU dynamics are in delta
        % formulation or not
        function obj = compute_Ref_Constraints(obj, delta_config)
            [obj.di_ref, obj.Iti_ref] = compute_ref(obj.nb_subsystems, obj.Agraph,...
                                                    obj.Vr, obj.Il, obj.Rij_mat, obj.Ri, obj.Vin); 
            obj = obj.setConstraints(obj, delta_config);
           
        end
        
    end % end methods

    methods (Static)       
        function obj = setSelectionMatrices(obj)
               state_i = eye(obj.ni);
               N = obj.nb_subsystems*obj.ni;
               for i=1:obj.nb_subsystems
                    obj.U{i} = [zeros(obj.ni,obj.ni*(i-1)), state_i, zeros(obj.ni,N-obj.ni*i)];
                    out_neighbors = sort([i;neighbors(obj.NetGraph, i)]); % neighbor states of i
                    obj.W{i} = zeros(obj.nb_subsystems);
                    % put a 1 in diagonal of neighbor state
                    for j=1:length(out_neighbors)
                        obj.W{i}(out_neighbors(j),out_neighbors(j)) = 1;
                    end
                    obj.W{i}= kron(obj.W{i}, state_i); % each subsystem has size mi
                    obj.W{i}(all(obj.W{i}==0,2),:)=[]; % remove all zero rows
                end
                % Create Wij: extract state j from neighbor set of state i (j belongs to
                % xNi)
                for i=1:obj.nb_subsystems
                    out_neighbors = sort([i;neighbors(obj.NetGraph, i)]); % neighbor states of i
                    for j = 1:length(out_neighbors)
                         obj.Wij{i}{out_neighbors(j)} = obj.U{out_neighbors(j)}*obj.W{i}';
                    end
                end
        end
       
        % Set constraints, depending if the dynamics are in delta formulation or
        % not
        function obj = setConstraints(obj, delta_config)
          obj.delta_config = delta_config;
          for i= 1:obj.nb_subsystems
            obj.Gx_i{i}=  [eye(obj.ni); -eye(obj.ni)];
            obj.Gu_i{i} = [1;-1];
            if ~obj.delta_config
                obj.Xref{i} = [obj.Vr(i); obj.Iti_ref(i)-obj.Il(i)];
                obj.fx_i{i} = [obj.Vmax(i); obj.Imax(i)-obj.Il(i);...
                                -obj.Vmin(i); -obj.Imin(i)+obj.Il(i)];
                obj.Uref{i} = obj.di_ref(i)-obj.Ri(i)*obj.Il(i)/obj.Vin(i);
                obj.fu_i{i} = [1;0];            
            elseif obj.delta_config
                obj.Xref{i} = [obj.Vr(i); obj.Iti_ref(i)];
                obj.Uref{i} = obj.di_ref(i);
                obj.fx_i{i} = [obj.Vmax(i); obj.Imax(i); -obj.Vmin(i); ...
                          -obj.Imin(i)] + [-obj.Xref{i}; obj.Xref{i}];
                obj.fu_i{i} = [1;0]+[-obj.Uref{i}; obj.Uref{i}];

            end
          end 
          for i= 1:obj.nb_subsystems % neighbors constraints only for active DGUs
            out_neighbors = sort([i;neighbors(obj.NetGraph, i)]); % neighbor states of
            obj.Gx_Ni{i} = blkdiag(obj.Gx_i{out_neighbors});
            obj.fx_Ni{i} = vertcat(obj.fx_i{out_neighbors});
          end   
        end    
        
        % Plots the DGU states (Voltage & current) and the control input (duty
        % cycle)
        function plot_DGU_system(X,U, config, control_type, param, simStart,...
                            dgu2plot, annot2plot)  
         M = length(dgu2plot);
         lgd = cell(1,M);
         voltage = cell(1,M); current= cell(1,M); integrator = cell(1,M);
         if config == "GENERAL" % not distributed: plot states for general sys
            states = cell2mat(X); % extract position/velocity at each timesteps
            k = M * param.ni;
            j = 1;
            for i = dgu2plot
                voltage{i} = states(j:k:end);
                current{i} = states(j+1:k:end);
                integrator{i} = states(j+2:k:end);
                j = j+param.ni;
                lgd{i} = sprintf("DGU %d", i);
            end
            controller = cell2mat(U');%first row u1, second row u2, column are timesteps
        elseif config == "DISTRIBUTED"
            k = param.ni;
            states = cell2mat(X');
            controller = cell2mat(U')'; %first row u1, second row u2, column are timesteps
            for i = dgu2plot
                voltage{i} = states(1:k:end,i);
                current{i} = states(2:k:end,i);
                if ~param.delta_config
                    current{i} = current{i}+ repmat(param.Il(i), size(current{i}));
                    controller(i,:) = controller(i,:)+repmat(param.Il(i)*param.Ri(i)/...
                                                         param.Vin(i), size(controller(i,:)));
                end
                lgd{i} = sprintf("DGU %d", i);
            end
            
        else
            error("not implemented configuration in plot states");
         end
 
        sim_stepsX = (simStart-1):length(X)-1;
        tx = param.Ts .* sim_stepsX; % x-time vector for states in seconds
        sim_stepsU = (simStart-1):length(U)-1;
        tu = param.Ts .* sim_stepsU; % x-time vector for controller in seconds
        figure()
        sgtitle(control_type);
        subplot(2,1,1)
        % --VOLTAGES--%
        title('Converter Voltages');
        hold on
        for i = dgu2plot
            plot(tx,voltage{i}, 'linewidth', 1.2);
        end
        if exist('annot2plot') && ~isempty(annot2plot)
            xline(annot2plot.array, '--');
            text(annot2plot.array, [50.15 50.15], annot2plot.text);
        end
        legend(string(lgd(dgu2plot)), 'FontSize', 6);
        grid on
        ylabel('Voltage [V]');
        xlabel('Time [s]');
        hold off
        subplot(2,1,2)
        %--CURRENTS--%
        title('Converter Currents');
        hold on
        for i = dgu2plot
            plot(tx, current{i},'linewidth', 1.2);
        end
        if exist('annot2plot') && ~isempty(annot2plot)
             xline(annot2plot.array, '--');
            text(annot2plot.array, [2.5 2.5], annot2plot.text);
        end
        legend(string(lgd(dgu2plot)), 'FontSize', 6);
        ylabel('Current [A]');
        grid on
        xlabel('Time [s]');
        hold off

        figure()
        %--DUTY CYCLE--%
        title("Controller")
        hold on
        for i = dgu2plot
            plot(tu, controller(i,:),'linewidth', 1.2);
        end
        if exist('annot2plot') && ~isempty(annot2plot)
            xline(annot2plot.array, '--');
            text(annot2plot.array, [0.27 0.27], annot2plot.text);
            xlim([0;1e-4]);
        end
        legend(string(lgd(dgu2plot)), 'FontSize', 6);
        xlabel('Time [s]');
        ylabel('Duty cycle');
        grid on
        hold off
        end
       
    end  % end static methods
end  % end class

% Function to change from the laplacian graph based representation (used in the
% passivity based redesign phase) to the system xi+ = A_Ni x_Ni + Bi ui
function A_Ni = change_system_representation(Ai,Fi,Ci,Agraph)
            M = size(Ai,2); % number of subsystems    
            G = graph(Agraph);
            A_Ni = cell(1,M);
            for i=1:M
                out_neighbors = sort([i;neighbors(G, i)]);
                Acell = cell(1,length(out_neighbors));
                for j=1:length(out_neighbors)
                    if isequal(out_neighbors(j), i)
                        Acell{j} = Ai{i} - Fi{i}*sum(Agraph(i,:))*Ci{i};    
                    else
                        Acell{j} = Fi{i}*Agraph(i, out_neighbors(j))*...
                                             Ci{out_neighbors(j)};
                    end
                end
                A_Ni{i} = cell2mat(Acell);  
            end
end
% Computes input & current reference (setpoints/equilibrium) given 
% voltage references
function [di_ref, Iti_ref] = compute_ref(M, Agraph, Vr, Il, Rij, R, Vin)
            di_ref = zeros(1,M); Iti_ref = zeros(1,M);
            G = digraph(Agraph);
            for i= 1:M
                out_neighbors = sort([successors(G, i)]);
                Vi = repmat(Vr(i), size(Vr(out_neighbors)));
                sum_over_Ni = sum((Vi - Vr(out_neighbors))./Rij(i,out_neighbors));
                di_ref(i) = (Vr(i)+ Il(i)*R(i))/Vin(i) + R(i)/Vin(i)*sum_over_Ni;
                Iti_ref(i) = (di_ref(i)*Vin(i) - Vr(i))/R(i);
            end     
end
