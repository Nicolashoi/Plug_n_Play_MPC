

function A_Ni = change_system_representation(Ai,Fi,Ci,Agraph)
M = size(Ai,2); % number of subsystems    
G = digraph(Agraph);
for i=1:M
    out_neighbors = sort([i,successors(G, i)]);
    Acell = cell(1,length(out_neighbors));
    Acell{i} = Ai{i} - Fi{i}*sum(Agraph(i,:))*Ci{i};
    out_neighbors(i) = [];
    for j=1:length(out_neighbors)
       Acell{out_neighbors(j)} = Fi{i}*Agraph(i, out_neighbors(j))*...
                                 Ci{out_neighbors(j)};
    end
    A_Ni{i} = cell2mat(Acell);    
end