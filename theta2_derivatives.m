function [ddtheta2_dtheta2,alledges_bfs] = theta2_derivatives(pixelind_bfs, G, dualedgeInd, narrowBand)

pix_dualedgord = find(narrowBand==1);
alledges_bfs = sparse(size(dualedgeInd,1),2);
ddtheta2_dtheta2 = sparse(size(dualedgeInd,1), size(G.Nodes,1));

for nr_dualedge=1:size(dualedgeInd,1)
    alledges_bfs(nr_dualedge,1) = find(pix_dualedgord==dualedgeInd(nr_dualedge,1));
    alledges_bfs(nr_dualedge,2) = find(pix_dualedgord==dualedgeInd(nr_dualedge,2));
    ddtheta2_dtheta2(nr_dualedge, pixelind_bfs==alledges_bfs(nr_dualedge,2)) = 1;
    ddtheta2_dtheta2(nr_dualedge, pixelind_bfs==alledges_bfs(nr_dualedge,1)) = -1;
end


end

