
classdef utilityFunctions
    methods (Static)
          % Compute initial state and matrices 
          function [x0, Q_Ni, Ri, Qi] = tuningParam(dguNet, delta_config)
            Q_Ni = cell(1,dguNet.nb_subsystems); 
            Qi = cell(1, dguNet.nb_subsystems);
            Ri = cell(1,dguNet.nb_subsystems);
            x0 = cell(1,dguNet.nb_subsystems);
            for i = 1:dguNet.nb_subsystems
                if delta_config
                    x0{i} = [50;5]; % initial condition in normal coordinates
                    % delta formulation carried in the simulation functions
                    % directly
                elseif ~delta_config
                    x0{i} = [50;5-dguNet.Il(i)]; % second state is Ii - Il
                else
                    error("config delta must be true or false");
                end
                m_Ni = size(dguNet.W{i},1);
                Q_Ni{i} =1*eye(m_Ni);
                Ri{i} = 1*eye(size(dguNet.Bi{i},2));
                Qi{i} = eye(dguNet.ni);
            end
          end  
            
        % Import parameters for the DGU system from a .txt file  
        function [subsystems, Vin, R, L, C, Vmax, Vmin, Imax, Imin] = importData(filename)
            delimiterIn = ';';
            headerlinesIn= 1;
            dataDGU = importdata(filename, delimiterIn, headerlinesIn);
            subsystems = dataDGU.data(1,1);
            disp(dataDGU.colheaders);
            disp(dataDGU.data);
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
        function cost = compute_QR_cost(X,U,Q,R, config, Xref, Uref)
            cost = 0;
            xref = vertcat(Xref{:});
            uref = vertcat(Uref{:});
            for k=1:length(X)-1
                if config == "DISTRIBUTED"
                    x = reshape(X{k},[],1); u = U{k}';
                    cost = cost + (x-xref)'*Q*(x-xref) + ...
                                   (u-uref)'*R*(u-uref);
                elseif config == "GENERAL"
                    cost = cost + (X{k}-xref)'*Q*(X{k}-xref) + ...
                                   (U{k}-uref)'*R*(U{k}-uref);
                else 
                    disp("WARNING: wrong config, choose general or distributed");
                end
                    
            end
        end
        
        %% COmpute tracking error magnitude
        function err = tracking_error(X,U, config, param, dgu2compute)
            if config == "GENERAL" % not distributed: plot states for general sys
                states = cell2mat(X); % extract position/velocity at each timesteps
                k = M * param.ni;
                j = 1;
                for i = dgu2compute
                    voltage{i} = states(j:k:end);
                    current{i} = states(j+1:k:end);
                    integrator{i} = states(j+2:k:end);
                    j = j+param.ni;

                end
                controller = cell2mat(U');%first row u1, second row u2, column are timesteps
            elseif config == "DISTRIBUTED"
                k = param.ni;
                states = cell2mat(X');
                for i = dgu2compute
                    voltage{i} = states(1:k:end,i);
                    current{i} = states(2:k:end,i);
                    if ~param.delta_config
                        current{i} = current{i}+ repmat(param.Il(i), size(current{i}));
                    end
                end
                controller = cell2mat(U')'; %first row u1, second row u2, column are timesteps
            end
                
            % Compute tracking error here
            
        end
       
    end % end methods
end % end class