
classdef utilityFunctions
    methods (Static)
          % Compute initial state and matrices 
          function [x0, Q_Ni, Ri, Qi] = tuningParam(param, delta_config)
            Q_Ni = cell(1, length(param.activeDGU)); 
            Qi = cell(1, length(param.activeDGU));
            Ri = cell(1,length(param.activeDGU));
            x0 = cell(1,param.nb_subsystems);
            for i = 1:param.nb_subsystems
                if delta_config
                    x0{i} = [50;5]; % initial condition in normal coordinates
                    % delta formulation carried in the simulation functions
                    % directly
                elseif ~delta_config
                    x0{i} = [50;5-param.Il(i)]; % second state is Ii - Il
                else
                    error("config delta must be true or false");
                end
            end
            for i=param.activeDGU
                [Qi{i}, Ri{i}, decVariables{i}] = computeQi_Ri(param, i);
            end
            for i=param.activeDGU
                neighbors_i = sort([i; neighbors(param.NetGraph, i)]);
                Q_Ni{i} = blkdiag(Qi{neighbors_i});
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
        
        %% Compute tracking error magnitude
        function err = tracking_error(X,U, config, param, dgu2compute)
            % Compute tracking error here
            err = 0;
            for k = 1:length(X)-1
                if config == "DISTRIBUTED"
                    err = err + sum((X{k}(:, dgu2compute)- ...
                                horzcat(param.Xref{dgu2compute})).^2, 'all')...
                          + sum((U{k}(dgu2compute)'- ...
                                 vertcat(param.Uref{dgu2compute})).^2);
                elseif config == "GENERAL"
                    err = err + sum((X{k}- vertcat(param.Xref{:})).^2, 'all')...
                            + sum((U{k}-vertcat(param.Uref{:})).^2);
                end
            end
            err = sqrt(err);
        end
       
    end % end methods
end % end class