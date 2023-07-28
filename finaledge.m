function [H, dualedgeInd]=finaledge(M, BoundaryPixelInd, intEdgeInd, narrowBandEdges, numBoundaries)

%Each primal edge (taken from intEdgeInd) is rotated by 90 degrees in the
%counterclockwise direction to get its corresponding dual edge.
dualedgeInd=narrowBandEdges(intEdgeInd,:);
for i=1:length(dualedgeInd)
    Delta=dualedgeInd(i,2)-dualedgeInd(i,1);
    j = dualedgeInd(i,1);
    if Delta==M
        dualedgeInd(i,2)=j-1; 
    else 
        dualedgeInd(i,1)=j-1;
        dualedgeInd(i,2)=j-M-1;
    end
end

E=size(dualedgeInd,1);
H=zeros(numBoundaries,E);

X = cell(numBoundaries, 1);
for p=1:numBoundaries
   %Connecting consecutive boundary pixels to get the boundary dual edges
   OrientedBDE = [BoundaryPixelInd{p}, circshift(BoundaryPixelInd{p}, -1)];
   %Assigning 1 for boundary dual edges in the direction of the exterior (counterclockwise) and
   %interior boundaries (clockwise) and -1 for those against 
   X = ismember(dualedgeInd, [OrientedBDE(:,1), OrientedBDE(:,2)], 'rows');
   Y = ismember(dualedgeInd, [OrientedBDE(:,2), OrientedBDE(:,1)], 'rows');
   pHrow = X - Y;
   H(p,:) = pHrow;
end

%H = -H; % why?
