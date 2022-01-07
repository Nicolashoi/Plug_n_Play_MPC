
classdef utilityFunctions
    methods (Static)
          % Compute initial state and matrices 
          function [x0, Q_Ni, Ri, Qi, decVariables] = tuningParam(param, delta_config, passivity)
            M = param.nb_subsystems;
            Q_Ni = cell(1, M); 
            Qi = cell(1, M);
            Ri = cell(1,M);
            x0 = cell(1,M);
            %% Set Initial State
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
            %% Get Qi, Ri and Q_Ni
            % If passivity is used, it is necessary to find Qi and Ri satisfying
            % the global Riccati equation to ensure asymptotic stability
            switch passivity
                case true
                    for i=param.activeDGU
                        [Qi{i}, Ri{i}, decVariables{i}] = computeQi_Ri(param, i);
                    end
                    decVarSum = sum(cell2mat(decVariables),2);
                    if all(decVarSum <= 1)
                        fprintf(['Qi and Ri found to guarantee asympt. stability '...
                                 'of the global system \n']);
                    else
                        fprintf(['WARNING: no Qi and Ri found such that asympt. stability '...
                                 'of the global system is guaranteed \n']);
                    end
                case false
                        Qi(:) = {2*eye(param.ni)};
                        Ri(:) = {eye(param.nu)};
       
                otherwise 
                    disp('specify if passivity is used or not for redesign phase');
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
                    xref = horzcat(param.Xref{dgu2compute});
                    err = err + sum((X{k}(2, dgu2compute)- xref(2,dgu2compute)).^2, 'all');%...
%                           + sum((U{k}(dgu2compute)'- ...
%                                  vertcat(param.Uref{dgu2compute})).^2);
                elseif config == "GENERAL"
                    err = err + sum((X{k}- vertcat(param.Xref{:})).^2, 'all')...
                            + sum((U{k}-vertcat(param.Uref{:})).^2);
                end
            end
            err = sqrt(err);
        end
       
    end % end methods
end % end class