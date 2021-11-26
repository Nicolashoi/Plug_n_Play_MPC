%% Main File for PnP

%% Initial DGU Network
nb_subsystems = 6;
Rij_mat = zeros(nb_subsystems);
Rij_mat(1,2) = 1.75; Rij_mat(2,3) = 3.5; Rij_mat(2,4) = 1.75; 
Rij_mat(3,5) = 2.8; Rij_mat(3,6) = 1.75;
Rij_mat = Rij_mat + tril(Rij_mat',1);