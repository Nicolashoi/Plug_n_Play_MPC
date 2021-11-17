%% DGU Parameter file

R{1} = 3e-3; Li{1} = 92e-6; Ct{1} = 75e-6;
R{2} = 1.5e-3; Li{2} = 102e-6; Ct{2} = 50e-6;
%Rij{1} = 1.75; Rij{2} = Rij{1};
R12 = 1.75;
Rij = [0 R12; R12 0];
Agraph = [0 1/R12; 1/R12 0];
Vin{1} = 100; Vin{2} = 100;
Vr{1} = 50; Il{1} = 5; 
Vr{2} = 50; Il{2} = 5;
nb_subsystems = 2; 
L = diag(sum(Agraph))-Agraph;