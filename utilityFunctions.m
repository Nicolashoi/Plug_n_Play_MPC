
classdef utilityFunctions
    methods (Static)
        

        
      function [x0, Q_Ni, Ri] = tuningParam(dguNet, delta_config)
            Q_Ni = cell(1,dguNet.nb_subsystems); 
            Ri = cell(1,dguNet.nb_subsystems);
            x0 = cell(1,dguNet.nb_subsystems);
            for i = 1:dguNet.nb_subsystems
                if delta_config
                    x0{i} = [50;5];
                elseif ~delta_config
                    x0{i} = [50;0]; % second state is Ii - Il
                else
                    error("config delta must be true or false");
                end
                m_Ni = size(dguNet.W{i},1);
                Q_Ni{i} =1*eye(m_Ni);
                Ri{i} = 1*eye(size(dguNet.Bi{i},2));
            end
        end

        function [subsystems, Vin, R, L, C, Vmax, Vmin, Imax, Imin] = importData(filename)
            delimiterIn = ';';
            headerlinesIn= 1;
            dataDGU = importdata(filename, delimiterIn, headerlinesIn);
            subsystems = dataDGU.data(1,1);
            Vin = dataDGU.data(:,2);
            R = dataDGU.data(:,3);
            L = dataDGU.data(:,4);
            C = dataDGU.data(:,5);
            Vmax = dataDGU.data(:,6);
            Vmin = dataDGU.data(:,7);
            Imax = dataDGU.data(:,8);
            Imin = dataDGU.data(:,9);
        end 
       
        %% Compute Finite Cost given state, input, Q and R
        function cost = compute_QR_cost(X,U,Q,R, config)
            cost = 0;
            if config == "GENERAL"
                for i=1:length(X)-1
                    cost = cost + X{i}'*Q*X{i} + U{i}'*R*U{i};
                end
            elseif config == "DISTRIBUTED"
                for i=1:length(X)-1
                    x = reshape(X{i},[],1); u = U{i}';
                    cost = cost + x'*Q*x + u'*R*u;
                end
            else
                error("wrong configuration, choose between general or distributed");
            end   
        end
       
    end % end methods
end % end class