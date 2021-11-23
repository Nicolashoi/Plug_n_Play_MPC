

classdef DGU_network
    properties
       %% Electrical Parameters
       nb_subsystems
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
        
    end
    
    methods
        function obj = DGU_network(nb_subsystems)
            if nargin == 1
                obj.nb_subsystems = nb_subsystems;
            end
        end
        function obj = initElecParam(obj, idx_DGU, Vin, Vr, Il, Ri, Ct, Li)
            obj.Vin(idx_DGU) = Vin;
            obj.Vr(idx_DGU) = Vr;
            obj.Il(idx_DGU) = Il;
            obj.Ri(idx_DGU) = Ri;
            obj.Ct(idx_DGU) = Ct;
            obj.Li(idx_DGU) = Li;
        end
            
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
                % compute references here
            end
            
           
        end
        
    end
    
    
    
end