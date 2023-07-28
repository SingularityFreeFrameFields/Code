function [theta1_BFS_order, theta1_matrix]=integrate_theta(phi, d_theta1, edges_bfs, roots_bfs, pixelind_bfs, Tree, G, dualedgeInd, narrowBand)

numPixelsNarrowBand = size(G.Nodes,1);
target = edges_bfs(:,2);
source = edges_bfs(:,1);
sourcecol= G.Nodes.PixelIndex(source)';
targetcol= G.Nodes.PixelIndex(target)';
finalgedge = [sourcecol;targetcol]';

[~,ia1,ib1] = intersect(finalgedge, dualedgeInd, "rows");
flipdualedgeInd=[dualedgeInd(:,2)';dualedgeInd(:,1)']';
[~,ia2,ib2] = intersect(finalgedge, flipdualedgeInd, "rows");
signs = -ones(size(finalgedge,1),1); %signs of bfs edges
signs(ia1) = 1;

dualEdges_correctOrder = zeros(size(finalgedge,1),1);
dualEdges_correctOrder(ia1) = ib1;
dualEdges_correctOrder(ia2) = ib2;
angles = d_theta1(dualEdges_correctOrder).*signs;

theta1_graphVertexOrder = zeros(numPixelsNarrowBand,1);
theta1_graphVertexOrder(roots_bfs) = phi;

for i = 1:size(edges_bfs,1)
    s = edges_bfs(i,1);
    t = edges_bfs(i,2);
    theta1_graphVertexOrder(t) = theta1_graphVertexOrder(s)+ angles(i);
end
theta1_BFS_order = theta1_graphVertexOrder(pixelind_bfs);

finalinfo_theta1 = table2array(Tree.Nodes);
angleforvec1 = zeros(size(narrowBand));
thirdcolumn1 = finalinfo_theta1(:,3);

for theta1_index=1:size(pixelind_bfs,1)-size(roots_bfs,2)
    angleforvec1(thirdcolumn1(pixelind_bfs(theta1_index)))=theta1_BFS_order(theta1_index);
end
theta1_matrix = angleforvec1;