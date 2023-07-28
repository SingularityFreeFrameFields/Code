function [dtheta1_ddtheta1,ddtheta2_ddtheta1] = theta1_derivatives(roots_bfs,edges_bfs, pixelind_bfs, G, dualedgeInd)
dtheta1_ddtheta1= sparse(size(G.Nodes,1),size(dualedgeInd,1));

target = edges_bfs(:,2);
source = edges_bfs(:,1);
dtheta1_ddtheta1(roots_bfs,:)=0;

sourcecol= G.Nodes.PixelIndex(source)';
targetcol= G.Nodes.PixelIndex(target)';
finalgedge = [sourcecol;targetcol]';

[~,ia1,ib1] = intersect(finalgedge, dualedgeInd, "rows");
flipdualedgeInd=[dualedgeInd(:,2)';dualedgeInd(:,1)']';
[~,ia2,ib2] = intersect(finalgedge, flipdualedgeInd, "rows");
dualEdgeIndex_toBFSIndex = zeros(size(dualedgeInd,1),1);
dualEdgeIndex_toBFSIndex(ib1) = ia1;
dualEdgeIndex_toBFSIndex(ib2) = ia2;
c_oneind_dtheta1_ddtheta1 = [];

%first assume that the columns of the matrix correspond to edge indices in
%BFS order. reorder after the loop
testIndPos = cell(size(G.Nodes,1),1);
testIndNeg = cell(size(G.Nodes,1),1);

iInd = [size(G.Nodes,1)];
jInd = [size(dualedgeInd,1)];
val = [0];
for i=1:size(edges_bfs,1)
 s = source(i);
 t = target(i);

 if (ismember(i,ia1))
     testIndPos{t} = [testIndPos{s};i];
     testIndNeg{t} = testIndNeg{s};
 else
     testIndPos{t} = testIndPos{s};
     testIndNeg{t} = [testIndNeg{s};i];
 end
 iInd = [iInd; ones(numel(testIndPos{t})+numel(testIndNeg{t}),1)*t];
 jInd = [jInd; testIndPos{t}; testIndNeg{t}];
 val = [val; ones(numel(testIndPos{t}),1); -ones(numel(testIndNeg{t}),1)];
end

dtheta1_ddtheta1 = sparse(iInd, jInd, val);

edges_not_in_BFS = find(~dualEdgeIndex_toBFSIndex);
numEdgesInBFS = size(edges_bfs,1);
numEdgesTotal = size(dualedgeInd,1);
dualEdgeIndex_toBFSIndex(edges_not_in_BFS) = [numEdgesInBFS+1:numEdgesTotal]'; %for the edges that are not in BFS, fill in dummy indices
dtheta1_ddtheta1 = dtheta1_ddtheta1(:, dualEdgeIndex_toBFSIndex);
dtheta1_ddtheta1 = dtheta1_ddtheta1(pixelind_bfs,:); %convert to BFS order of rows

ddtheta2_ddtheta1 = sparse(size(dualedgeInd,1), size(dualedgeInd,1));

for i=1:size(dualedgeInd,1)
    ddtheta2_ddtheta1(i,i)=1;
end

end

