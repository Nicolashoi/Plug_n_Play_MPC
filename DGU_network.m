

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
        
       %% Others
       Agraph
       graph
       Rij_mat
       Iti_ref
       di_ref
       nb_subsystems
       ni = 2 % size of subsystem state
       %% Constraints
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
    end
    
    methods
        % constructor
        function obj = DGU_network(nb_subsystems, Rij_mat)
            if nargin == 2
                obj.nb_subsystems = nb_subsystems;
                obj.Rij_mat = Rij_mat;
                obj.Agraph = Rij_mat;
                nonzeroIdx = Rij ~= 0;
                obj.Agraph(nonzeroIdx) = 1./Agraph(nonzeroIdx);  
                obj.graph = digraph(obj.Agraph);
            end
        end
        
        % initialize electrical parameters
        function obj = initElecParam(obj, idx_DGU, Vin, Vr, Il, Ri, Ct, Li)
            obj.Vin(idx_DGU) = Vin;
            obj.Vr(idx_DGU) = Vr;
            obj.Il(idx_DGU) = Il;
            obj.Ri(idx_DGU) = Ri;
            obj.Ct(idx_DGU) = Ct;
            obj.Li(idx_DGU) = Li;
        end
        
        % initialize DGU dynamics
        function obj = initDynamics(obj,idx_DGU)
            for i = idx_DGU  % index of DGU (can be multiple) to get dynamics
                if isempty(obj.Ct) || isempty(obj.Vr)
                    error("Electrical parameters not correctly initialized or missing")
                end
                Ac = [0, 1/obj.Ct(i); -1/obj.Li(i), -obj.Ri(i)/obj.Li(i)];
                Bc = [0; obj.Vin(i)/obj.Li(i)];
                Fc = [1/obj.Ct(i); 0];
                Cc = [1, 0];  
                subsys = ss(Ac, [Bc Fc], [], []);
                subsys_d = c2d(subsys, obj.Ts);
                obj.Ai{i} = subsys_d.A;
                obj.Bi{i} = subsys_d.B(:,1);
                obj.Fi{i} = subsys_d.B(:,2);
                obj.Ci{i} = Cc;     
            end
            % function to change representation to A_Ni
            obj.A_Ni = change_system_representation(obj.Ai, obj.Fi, obj.Ci, obj.Agraph);
            % function compute references here
            [obj.di_ref, obj.Iti_ref] = compute_ref(obj.nb_subsystems, obj.Agraph,...
                                                    obj.Vr, obj.Il, obj.Rij, obj.Ri, obj.Vin)    
            % Compute constraints
            set_constraints()
        end
        
        % change representation to neighbor state matrix
        function obj = change_system_representation(obj)   
            %A_Ni = cell(1,M);
            if length(obj.Ai) ~= obj.nb_subsystems
                error("Provide all system dynamics before changing to A_Ni");
            end
            for i=1:obj.nb_subsystems
                out_neighbors = sort([i;successors(obj.graph, i)]);
                Acell = cell(1,length(out_neighbors));
                for j=1:length(out_neighbors)
                    if isequal(out_neighbors(j), i)
                        Acell{j} = obj.Ai{i} - obj.Fi{i}*sum(obj.Agraph(i,:))*obj.Ci{i};
                    else
                        Acell{j} = obj.Fi{i}*obj.Agraph(i, out_neighbors(j))*...
                                             obj.Ci{out_neighbors(j)};
                    end
                end
                obj.A_Ni{i} = cell2mat(Acell);  
            end
        end
         
        function obj = compute_ref(obj)
            %di_ref = zeros(1,M); Iti_ref = zeros(1,M);
            G = digraph(obj.Agraph);
            for i= 1:obj.nb_subsystems
                out_neighbors = sort(successors(G, i));
                Vi = repmat(obj.Vr(i), size(obj.Vr(out_neighbors)));
                sum_over_Ni = sum((Vi - obj.Vr(out_neighbors))./...
                              obj.Rij_mat(i,out_neighbors));
                obj.di_ref(i) = (obj.Vr(i)+ obj.Il(i)*obj.Ri(i))/obj.Vin(i) + ...
                                 obj.Ri(i)/obj.Vin(i)*sum_over_Ni;
                obj.Iti_ref(i) = (obj.di_ref(i)*obj.Vin(i) - obj.Vr(i))/obj.Ri(i);
            end     
        end
        
        function obj = setConstraints(obj)
              for i= 1:obj.nb_subsystems
                obj.Xref{i} = [obj.Vr(i); obj.Iti_ref(i)-obj.Il(i)];
                obj.Uref{i} = obj.di_ref{i};
                obj.Gx_i{i}=  [eye(2); -eye(2)];
                obj.Gu_i{i} = [1;-1];
                obj.fx_i{i} = [52; 10-obj.Il(i); -49; 0+obj.Il(i)];
                obj.fu_i{i} = [1;0];
              end 
              for i= 1:obj.nb_subsystems
                out_neighbors = sort([i;successors(obj.graph, i)]); % neighbor states of
                obj.Gx_Ni{i} = blkdiag(obj.Gx_i{out_neighbors});
                obj.fx_Ni{i} = vertcat(obj.fx_i{out_neighbors});
              end   
        end
        
    end  % end methods
end  % end class


function A_Ni = change_system_representation(Ai,Fi,Ci,Agraph)
            M = size(Ai,2); % number of subsystems    
            G = digraph(Agraph);
            A_Ni = cell(1,M);
            for i=1:M
                out_neighbors = sort([i;successors(G, i)]);
                Acell = cell(1,length(out_neighbors));
                %Acell{i} = Ai{i} - Fi{i}*sum(Agraph(i,:))*Ci{i};
                %out_neighbors(i) = [];
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

function [di_ref, Iti_ref] = compute_ref(M, Agraph, Vr, Il, Rij, R, Vin)
            %load('system/DGU_electrical_param.mat')
            %M = nb_subsystems;
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
        