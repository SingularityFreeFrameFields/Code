function [OldIndices,d_theta,singsxy,singIndex,path_per_signularity,backwardsLinearNarrowBand,preclean_dtheta] = CleanInitialization(Init, narrowBand, d0sub, H, intEdges, primalGraph, NBprimalGraph, intVertices, numBoundaries,...
    numExteriorBoundaries, smallLenBoundaryInd, BoundaryPixelInd, numExteriorBoundariesOld)
%% Docstring
%{
  Takes in a field initialization obtained through a simple convex
  optimization and returns a set of indices per boundary that may lead to
  nicer frame fields.
%}

%% Variable Setup

[M, N] = size(narrowBand);

% Complex Number representation of vector field
ComplexInit = Init(1:length(Init) / 2) + 1i * Init(length(Init) / 2 + 1 : end);

%{ 
  Construct a Graph from the NarrowBand in this code pixels are considered 
  DUAL VERTICES Nodes == Black Pixels Edges == Connections between Black Pixels
%}
G = binaryImageGraph(narrowBand, 4);

% Linear Indices of black pixels of NarrowBand
linearNarrowBand = G.Nodes.PixelIndex;

%{
 backwardsLinearNarrowBand(linearNarrowBand) = 1:numel(narrowBand) In
 other words, returns the implicit ordering of a narrowBand pixel Ex:
 backwardsLinearNarrowBand(50) == 16 This means the pixel with linear
 index '50' is the 16th pixel in the narrowBand
%}
backwardsLinearNarrowBand = zeros(numel(narrowBand),1);
backwardsLinearNarrowBand(linearNarrowBand) = (1:length(linearNarrowBand))';

%% Uncomment to see the Complex Number Representation of Field Initialization
%{
[i,j] = find(narrowBand);
quiver(j, -i, real(ComplexInit), -imag(ComplexInit), 0.1, 'red');
%}
%% Calculate d(theta) for NarrowBand

% Find all vertical and horizontal Primal Edges
vertPrimalEdges = find(intEdges(:,1) - intEdges(:,2) == 1);
horPrimalEdges = find(intEdges(:,1) - intEdges(:,2) ~= 1);

% save interior edges indices and ordering for later
intEdgesOrig = intEdges;

% Flip primal edges so that they correspond to dual edges
% OPPOSITE CONVENTION; CANCELLED A FEW LINES LATER
intEdges(vertPrimalEdges,1) = intEdges(vertPrimalEdges, 1) - M - 1;
intEdges(horPrimalEdges, 1) = intEdges(horPrimalEdges, 1) - 1;
intEdges(horPrimalEdges, 2) = intEdges(horPrimalEdges, 2) - M;

PixelEdges = backwardsLinearNarrowBand(intEdges);
complex_edges = ComplexInit(PixelEdges);
d_theta = angle(complex_edges(:, 1) ./ complex_edges(:, 2)); % CANCELLED HERE

%% Finding Singulaties

%{ 
  When d0dtheta is non-zero, then there is a singularity at the primal
  vertex
%} 
d0dtheta = (d0sub)' * d_theta;
sings = find(abs(d0dtheta) > 1e-15);
singIndex = d0dtheta(sings) ./ (2 * pi); %maybe we want to permute sing index but not sure rn (as of 4/29) 
%[singularPrimalVerticesInd, ~] = find(d0sub(:, sings));

% Translate the indexing of primal vertices used in buildDerivatives
% to the "real" primal graph indexing
singularityV = intVertices(sings);
[rowSV, colSV] = ind2sub([M, N], singularityV);
singularities = singularityV + (colSV - 1);

posSingsLog = singIndex > 0;
negSingsLog = singIndex < 0;

posSingsToSings = find(posSingsLog);
negSingsToSings = find(negSingsLog);

singsPlus = singularities(posSingsLog);
singsMinus = singularities(negSingsLog);

singsxy = [rowSV, colSV];

% Just a placeholder to make MATLAB happy that this is being assigned
path_per_signularity = [];

%% Creating relevant digraphs for transporting plus and minus singularities

% Create a digraph
NBprimalDigraph = digraph(adjacency(NBprimalGraph));
NBplusPrimalDigraph = NBprimalDigraph;
NBminusPrimalDigraph = NBprimalDigraph;

% Store edge endpoints and weights copies
primalDigraphEdges = NBprimalDigraph.Edges.EndNodes;
% Ensure boundary edges are more expensive than interior edges
newWtsPlus = 10*pi*NBprimalDigraph.Edges.Weight;
newWtsMinus = 10*pi*NBprimalDigraph.Edges.Weight;

% Convert to primal indices into top-left corner grid
[~,colPDE] = ind2sub([M+1,N+1],primalDigraphEdges);
TLprimalDigraphEdges = primalDigraphEdges - (colPDE-1);

% Find edges that agree with standard orientation
[withInds,withLocs] = ismember(TLprimalDigraphEdges,intEdgesOrig,'rows');
withLocs = nonzeros(withLocs);
% Assign the change in smoothness (up to constant scaling)
% newWtsPlus(withInds) = pi - d_theta(withLocs);
% newWtsMinus(withInds) = pi + d_theta(withLocs);
% Trying to incorporate initialization norm as well
newWtsPlus(withInds) = (abs(complex_edges(withLocs,1)) + abs(complex_edges(withLocs,2))).*(pi - d_theta(withLocs));
newWtsMinus(withInds) = (abs(complex_edges(withLocs,1)) + abs(complex_edges(withLocs,2))).*(pi + d_theta(withLocs));
% Do same for edges that disagree
[againstInds,againstLocs] = ismember(circshift(TLprimalDigraphEdges,1,2),intEdgesOrig,'rows');
againstLocs = nonzeros(againstLocs);
% Assign the change in smoothness (up to constant scaling)
% newWtsPlus(againstInds) = pi + d_theta(againstLocs);
% newWtsMinus(againstInds) = pi - d_theta(againstLocs);
% Trying to incorporate initialization norm as well
newWtsPlus(againstInds) = (abs(complex_edges(againstLocs,1)) + abs(complex_edges(againstLocs,2))).*(pi + d_theta(againstLocs));
newWtsMinus(againstInds) = (abs(complex_edges(againstLocs,1)) + abs(complex_edges(againstLocs,2))).*(pi - d_theta(againstLocs));

NBplusPrimalDigraph.Edges.Weight = newWtsPlus;
NBminusPrimalDigraph.Edges.Weight = newWtsMinus;

%% Computing Distances (old way)

% distancestoBoundary = [];
% closestVtx = [];
% 
% %{
% b = BoundaryPixelInd{1};
% [~, colB] = ind2sub([M, N], b);
% tb = BoundaryPixelInd{1} + (colB - 1);
% tB = [tb; tb + 1; tb + (M+1); tb + (M+2)];
% 
% [~, colIntV] = ind2sub([M,N], intVertices);
% tIntV = intVertices + (colIntV - 1);
% 
% tB = setdiff(tB, tIntV);
% %}
% 
% Uncomment to plot image with singularites labelled

% h = plot(primalGraph);
% %highlight(h, NBprimalGraph.Edges.EndNodes(:,1), NBprimalGraph.Edges.EndNodes(:,2), 'EdgeColor', 'black', 'NodeColor', 'red', 'LineWidth', 5);
% %highlight(h, tB, 'NodeColor', 'g');
% %highlight(h, tB, 'NodeColor', 'g', 'EdgeColor', 'g');
% highlight(h, singularities, 'NodeColor', 'm', 'MarkerSize', 5);
% axis equal
% set(gca, 'XDir','reverse')
% camroll(90);
% axis equal;
% axis off;
% axis equal off
% highlight(h, singularities, 'NodeColor', 'b', 'MarkerSize', 5);
% highlight(h, singularities, 'NodeColor', 'y', 'MarkerSize', 5);
% 
% 
% boundCoordsCell = getExactBounds(narrowBand);
% 
% boundCoordsCell_permuted = boundCoordsCell; 
% 
% numBoundariesold = length(BoundaryPixelIndold);
% 
% for i=1:numBoundariesold
%     firstPixel = min(BoundaryPixelIndold{i});
%     rows = size(narrowBand,1);
%     columns = size(narrowBand,2);
%     [firstPixelRow, firstPixelCol] = ind2sub([rows, columns], firstPixel);
%     for j=1:length(boundCoordsCell) %boundcoordscell outputs (x,y)-format coordinates; x is the col index, y is the row index 
%         if i <= numExteriorBoundaries
%             if ismember([firstPixelCol-0.5,firstPixelRow-0.5],boundCoordsCell{j},'rows')
%                 boundCoordsCell_permuted{i} = boundCoordsCell{j};
%                 break
%             end  
%         elseif ismember([firstPixelCol+0.5,firstPixelRow+0.5],boundCoordsCell{j},'rows')
%             boundCoordsCell_permuted{i} = boundCoordsCell{j};
%             break
%         end
%     end 
% end 
% 
% boundCoordsCell = boundCoordsCell_permuted; 
% 
% boundCoordsCell(smallLenBoundaryInd) = [];
% 
% primalBoundaryIdxCell = {size(boundCoordsCell)};
% 
% for i = 1:length(boundCoordsCell)
%     boundCoords = ceil(boundCoordsCell{i});
%     primalBoundaryIdxCell{i} = sub2ind([M+1, N+1], boundCoords(:,2), boundCoords(:,1));
% end
% 
% for i = 1:length(primalBoundaryIdxCell)
%     %[~, colB] = ind2sub([M, N], BoundaryPixelInd{i});
%     %translatedBoundaryPixels = BoundaryPixelInd{i} + (colB - 1);
%     %boundaryPixels = setdiff(translatedBoundaryPixels, tIntV);
%     %if (length(boundaryPixels) ~= 1)
%     primBound = unique(primalBoundaryIdxCell{i});
%     dist = distances(NBprimalGraph, singularities, primBound);
%     [distancestoBoundary(:, i),closestIdx] = min(dist, [], 2);
%     closestVtx(:,i) = primBound(closestIdx);
%     %end
% end
% 
% %[~, BoundariesToPush] = ind2sub(size(distancestoBoundary), I(M > 0));
% [~,I] = min(distancestoBoundary, [], 2);
% closestVtx = closestVtx(sub2ind(size(closestVtx),1:size(I),I'))';
% path_per_signularity = cell(size(singularities));
% for i=1:size(singularities)
%     [x,y] = ind2sub([M+1, N+1],shortestpath(NBprimalGraph,singularities(i),closestVtx(i)));
%     path_per_signularity{i} = [y',x'];
% end
% % %[~, BoundariesToPush] = ind2sub(size(distancestoBoundary), I(M > 0));

%% Computing Distances (Weighted Approach)

% Computing distances for moving minus singularities to plus singularities
minusToPlusDists = distances(NBminusPrimalDigraph, singsMinus, singsPlus);
% Returns x-by-0 or 0-by-x if one of the sets of singularities are empty

% Now for distances to the boundary
distancesMinusToBoundary = [];
closestMinusVtx = [];
distancesPlusToBoundary = [];
closestPlusVtx = [];

% A bunch of stuff for getting primal vertices as sets of boundaries
% boundCoordsCellOld = getExactBounds(narrowBand);

boundCoordsCell = {};

for i=1:numBoundaries
    [bdyLength,~] = size(BoundaryPixelInd{i});
    [bdyPixelRows,bdyPixelCols] = ind2sub([M,N],BoundaryPixelInd{i});
    bdyPixelXY = [bdyPixelRows,bdyPixelCols];
    bdyPixelXY = bdyPixelXY(1:(bdyLength-1),:); % get rid of repeated element
    % Tracing is different for exterior vs. interior boundaries
    if i <= numExteriorBoundaries
        if bdyLength == 2 % if the component is a single pixel
            topLeft = bdyPixelXY(1,:) - [0.5, 0.5];
            botLeft = bdyPixelXY(1,:) - [-0.5, 0.5];
            botRight = bdyPixelXY(1,:) - [-0.5, -0.5];
            topRight = bdyPixelXY(1,:) - [0.5, -0.5];
            boundCoordsCell{i,1} = [topLeft; botLeft; botRight; topRight];
        else
            [~,minInd] = min(BoundaryPixelInd{i});
            % Shift so minInd is at 1, for ease of conditions
            bdyPixelXY = circshift(bdyPixelXY,1-minInd, 1);
            boundCoordsCell{i,1} = bdyPixelXY(1,:)-[0.5, 0.5]; % initialize boundCoordsCell
            for j=2:bdyLength % recall first pixel repeated in output of bwboundaries
                if j == bdyLength
                    nextXY = bdyPixelXY(1,:);
                else
                    nextXY = bdyPixelXY(j,:);
                end
                lastXY = bdyPixelXY(j-1,:);
                if j==2
                    lastLastXY = bdyPixelXY(bdyLength-1,:);
                else
                    lastLastXY = bdyPixelXY(j-2,:);
                end
                prevStep = nextXY-lastXY;
                prevStep = mat2str(prevStep); % convert to char
                prevPrevStep = lastXY - lastLastXY;
                prevPrevStep = mat2str(prevPrevStep); % convert to char
                switch prevStep
                    case '[1 0]' % down to next pixel
                        switch prevPrevStep % note nothing happens for prev pixel from left
                            case '[0 -1]' % prev pixel from right
                                if j==2 % special start case
                                    newCorner = nextXY+[-0.5,-0.5]; 
                                    boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                                else
                                    newCorner1 = lastXY + [-0.5, -0.5];
                                    newCorner2 = nextXY + [-0.5, -0.5];
                                    boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                                end
                            case '[1 0]' % prev pixel from above
                                newCorner = nextXY+[-0.5,-0.5]; 
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                            case '[-1 0]' % prev pixel from below
                                if j==2 % special start case
                                    newCorner1 = lastXY + [-0.5, 0.5];
                                    newCorner3 = nextXY + [-0.5, -0.5];
                                    boundCoordsCell{i} = [newCorner1; boundCoordsCell{i}; newCorner3];
                                else
                                    newCorner1 = lastXY + [-0.5, 0.5];
                                    newCorner2 = lastXY + [-0.5, -0.5];
                                    newCorner3 = nextXY + [-0.5, -0.5];
                                    boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2; newCorner3];
                                end
                        end      
                    case '[0 1]' % right to next pixel
                        switch prevPrevStep % note nothing happens for prev pixel from below
                            case '[0 1]' % prev pixel from left
                                newCorner = nextXY + [0.5, -0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                            case '[1 0]' % prev pixel from above
                                newCorner1 = lastXY + [0.5, -0.5];
                                newCorner2 = nextXY + [0.5, -0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                            case '[0 -1]' % prev pixel from right
                                if j==2 % special start case
                                    newCorner1 = lastXY + [0.5, -0.5];
                                    newCorner2 = nextXY + [0.5, -0.5];
                                    boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                                else
                                    newCorner1 = lastXY + [-0.5, -0.5];
                                    newCorner2 = lastXY + [0.5, -0.5];
                                    newCorner3 = nextXY + [0.5, -0.5];
                                    boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2; newCorner3];
                                end
                        end     
                    case '[-1 0]' % next pixel above
                        switch prevPrevStep % note nothing happens for prev pixel from right
                            case '[-1 0]' % prev pixel from below
                                newCorner = nextXY + [0.5, 0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                            case '[0 1]' % prev pixel from left
                                newCorner1 = lastXY + [0.5, 0.5];
                                newCorner2 = nextXY + [0.5, 0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                            case '[1 0]' % prev pixel from above
                                newCorner1 = lastXY + [0.5, -0.5];
                                newCorner2 = lastXY + [0.5, 0.5];
                                newCorner3 = nextXY + [0.5, 0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2; newCorner3];
                        end     
                    case '[0 -1]' % next pixel to left
                        switch prevPrevStep % note nothing happens for prev pixel from above
                            case '[0 -1]' % prev pixel from right
                                newCorner = nextXY + [-0.5, 0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                            case '[-1 0]' % prev pixel from below
                                newCorner1 = lastXY + [-0.5, 0.5];
                                newCorner2 = nextXY + [-0.5, 0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                            case '[0 1]' % prev pixel from left
                                newCorner1 = lastXY + [0.5, 0.5];
                                newCorner2 = lastXY + [-0.5, 0.5];
                                newCorner3 = nextXY + [-0.5, 0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2; newCorner3];
                        end
                end
            end
        end
    else %interior boundary
        [~,minInd] = min(BoundaryPixelInd{i});
        % Shift so minInd is at 1, for ease of conditions
        bdyPixelXY = circshift(bdyPixelXY,1-minInd, 1);
        boundCoordsCell{i,1} = bdyPixelXY(1,:)+[0.5, 0.5]; % initialize boundCoordsCell
        for j=2:bdyLength % recall first pixel repeated in output of bwboundaries
            if j == bdyLength
                nextXY = bdyPixelXY(1,:);
            else
                nextXY = bdyPixelXY(j,:);
            end
            lastXY = bdyPixelXY(j-1,:);
            if j==2
                lastLastXY = bdyPixelXY(bdyLength-1,:);
            else
                lastLastXY = bdyPixelXY(j-2,:);
            end
            prevStep = nextXY-lastXY;
            prevStep = mat2str(prevStep); % convert to char
            prevPrevStep = lastXY - lastLastXY;
            prevPrevStep = mat2str(prevPrevStep); % convert to char
            switch prevStep
                case '[1 0]' % down to next pixel
                    switch prevPrevStep % note nothing happens for prev pixel from left
                        case '[0 -1]' % prev pixel from right
                            newCorner1 = lastXY + [-0.5, -0.5];
                            newCorner2 = nextXY + [-0.5, -0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                        case '[1 0]' % prev pixel from above
                            newCorner = nextXY+[-0.5,-0.5]; 
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                        case '[-1 0]' % prev pixel from below
                            newCorner1 = lastXY + [-0.5, 0.5];
                            newCorner2 = lastXY + [-0.5, -0.5];
                            newCorner3 = nextXY + [-0.5, -0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2; newCorner3];
                    end      
                case '[0 1]' % right to next pixel
                    switch prevPrevStep % note nothing happens for prev pixel from below
                        case '[0 1]' % prev pixel from left
                            if j ~= 2
                                newCorner = nextXY + [0.5, -0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                            end
                        case '[1 0]' % prev pixel from above
                            newCorner1 = lastXY + [0.5, -0.5];
                            newCorner2 = nextXY + [0.5, -0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                        case '[0 -1]' % prev pixel from right
                            newCorner1 = lastXY + [-0.5, -0.5];
                            newCorner2 = lastXY + [0.5, -0.5];
                            newCorner3 = nextXY + [0.5, -0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2; newCorner3];
                    end     
                case '[-1 0]' % next pixel above
                    switch prevPrevStep % note nothing happens for prev pixel from right
                        case '[-1 0]' % prev pixel from below
                            if j ~= bdyLength
                                newCorner = nextXY + [0.5, 0.5];
                                boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                            end
                        case '[0 1]' % prev pixel from left
                            newCorner1 = lastXY + [0.5, 0.5];
                            newCorner2 = nextXY + [0.5, 0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                        case '[1 0]' % prev pixel from above
                            newCorner1 = lastXY + [0.5, -0.5];
                            newCorner2 = lastXY + [0.5, 0.5];
                            newCorner3 = nextXY + [0.5, 0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2; newCorner3];
                    end     
                case '[0 -1]' % next pixel to left
                    switch prevPrevStep % note nothing happens for prev pixel from above
                        case '[0 -1]' % prev pixel from right
                            newCorner = nextXY + [-0.5, 0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner];
                        case '[-1 0]' % prev pixel from below
                            newCorner1 = lastXY + [-0.5, 0.5];
                            newCorner2 = nextXY + [-0.5, 0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2];
                        case '[0 1]' % prev pixel from left
                            newCorner1 = lastXY + [0.5, 0.5];
                            newCorner2 = lastXY + [-0.5, 0.5];
                            newCorner3 = nextXY + [-0.5, 0.5];
                            boundCoordsCell{i} = [boundCoordsCell{i}; newCorner1; newCorner2; newCorner3];
                    end
            end
        end
    end
    boundCoordsCell{i} = circshift(boundCoordsCell{i},1,2);
end

% Debugging visualization of the boundaries according to bwboundaries
% figure
% set(gca,'YDir','reverse');
% hold on;
% axis equal;
% for n=1:numel(BoundaryPixelInd)
%     [i,j] = ind2sub(size(narrowBand),BoundaryPixelInd{n});
%     lineObject = plot(j,i);
%     color = get(lineObject, 'Color');
%     centerIndex = round(numel(i/2));
%     centralPt = [i(centerIndex); j(centerIndex)];
%     text(centralPt(2), centralPt(1),strcat("#",int2str(n)),'Color',color);
% end

% Debugging visualization of the boundaries according to getExactBounds
% figure
% hold on
% imshow(narrowBand)
% hold on
% for k=1:length(boundCoordsCell)
%     hold on
%     plot(boundCoordsCell{k}(:,1),boundCoordsCell{k}(:,2),'s-g','MarkerSize',4,'Linewidth',1.5);
%     hold on
% end

% boundCoordsCell_permuted = boundCoordsCell; 

% TRYING TO GET RID OF SMALL LENGTH CAP
%numBoundariesold = length(BoundaryPixelIndold);

% for i=1:numBoundaries
%     firstPixel = min(BoundaryPixelInd{i});
%     rows = size(narrowBand,1);
%     columns = size(narrowBand,2);
%     [firstPixelRow, firstPixelCol] = ind2sub([rows, columns], firstPixel);
%     for j=1:length(boundCoordsCell) %boundcoordscell outputs (x,y)-format coordinates; x is the col index, y is the row index 
%         %if i <= numExteriorBoundariesOld
%         if i <= numBoundaries
%             if ismember([firstPixelCol-0.5,firstPixelRow-0.5],boundCoordsCell{j},'rows')
%                 boundCoordsCell_permuted{i} = boundCoordsCell{j};
%                 break
%             end  
%         elseif ismember([firstPixelCol+0.5,firstPixelRow+0.5],boundCoordsCell{j},'rows')
%             boundCoordsCell_permuted{i} = boundCoordsCell{j};
%             break
%         end
%     end 
% end 
% 
% boundCoordsCell = boundCoordsCell_permuted; 

%boundCoordsCell(smallLenBoundaryInd) = [];

primalBoundaryIdxCell = {size(boundCoordsCell)};

for i = 1:length(boundCoordsCell)
    boundCoords = ceil(boundCoordsCell{i});
    primalBoundaryIdxCell{i} = sub2ind([M+1, N+1], boundCoords(:,2), boundCoords(:,1));
end

% Let us do distances to the boundary for singularities
for i = 1:length(primalBoundaryIdxCell)
    %[~, colB] = ind2sub([M, N], BoundaryPixelInd{i});
    %translatedBoundaryPixels = BoundaryPixelInd{i} + (colB - 1);
    %boundaryPixels = setdiff(translatedBoundaryPixels, tIntV);
    %if (length(boundaryPixels) ~= 1)
    primBound = unique(primalBoundaryIdxCell{i});
    distMinus = distances(NBminusPrimalDigraph, singsMinus, primBound);
    distPlus = distances(NBplusPrimalDigraph, singsPlus, primBound);
    [distancesMinusToBoundary(:, i),closestMinusIdx] = min(distMinus, [], 2);
    closestMinusVtx(:,i) = primBound(closestMinusIdx);
    [distancesPlusToBoundary(:, i),closestPlusIdx] = min(distPlus, [], 2);
    closestPlusVtx(:,i) = primBound(closestPlusIdx);
    %end
end

[bdyDistsMinus,IMinus] = min(distancesMinusToBoundary, [], 2);
[bdyDistsPlus,IPlus] = min(distancesPlusToBoundary, [], 2);

% Assembling cost matrix
cost = [minusToPlusDists, bdyDistsMinus; bdyDistsPlus', 1e14];

% Solving assignment problem
costVectorized = cost(:);
costVectorized = costVectorized(1:(length(costVectorized)-1));
costVectorized(costVectorized == Inf) = 1e15; % linprog doesn't like Infs

% % Constructing row sums to 1 constraints
% Aeq = zeros(length(singsMinus)+length(singsPlus),length(costVectorized));
% for i=1:size(minusToPlusDists,1)
%     rowTemp = zeros(size(cost));
%     rowTemp(i,:) = 1;
%     rowTempVec = rowTemp(:);
%     rowTempVec = rowTempVec(1:(length(rowTempVec)-1));
%     Aeq(i,:) = rowTempVec';
% end
% % Constructing column sums to 1 constraints
% for i=1:size(minusToPlusDists,2)
%     colTemp = zeros(size(cost));
%     colTemp(:,i) = 1;
%     colTempVec = colTemp(:);
%     colTempVec = colTempVec(1:(length(colTempVec)-1));
%     Aeq(i+length(singsMinus),:) = colTempVec';
% end

% A sparse construction attempt
a = length(singsMinus);
b = length(singsPlus);
% Row sums
iRowSums = kron(1:a,ones(1,b+1));
linIndsCost = reshape(1:(a+1)*(b+1),[a+1,b+1]);
linIndsCostTrans = linIndsCost';
jRowSums = linIndsCostTrans(:)';
jRowSums = jRowSums(1:(length(costVectorized)-b));
% Column sums
iColSums = a+kron(1:b,ones(1,a+1));
jColSums = 1:(length(costVectorized)-a);
iSparse = [iRowSums, iColSums];
jSparse = [jRowSums, jColSums];
AeqSparse = sparse(iSparse,jSparse,ones(length(iSparse),1));
assn = linprog(costVectorized,[],[],AeqSparse,ones(size(AeqSparse,1),1),zeros(length(costVectorized),1),[]);

% put into matrix form for easy reading
assn = [assn; -1];
assn = reshape(assn,size(cost));
% closestVtx = closestVtx(sub2ind(size(closestVtx),1:size(I),I'))';
% path_per_signularity = cell(size(singularities));
% for i=1:size(singularities)
%     [x,y] = ind2sub([M+1, N+1],shortestpath(NBprimalGraph,singularities(i),closestVtx(i)));
%     path_per_signularity{i} = [y',x'];
% end
%% Create New Indices (old way)

% %ext 
% OldIndices(1: numExteriorBoundaries) = 1 - H(1: (numExteriorBoundaries), :)*d_theta ./ (2 * pi);
% %int 
% OldIndices((numExteriorBoundaries + 1): numBoundaries) = -1 - H((numExteriorBoundaries + 1): numBoundaries, :)*d_theta ./ (2 * pi);
% 
% for z=1:length(I)
%     OldIndices(I(z)) = OldIndices(I(z)) + singIndex(z);
% end
% 
% %{
% [count, bIndex] = groupcounts(BoundariesToPush);
% if (isempty(BoundariesToPush))
%    newIndices = OldIndices;
% else
%     OldIndices(bIndex) = OldIndices(bIndex) + count;
%     newIndices = OldIndices;
% end
% %}
% 
% %% Comb singularities out of dtheta (old way)
% 
% preclean_dtheta = d_theta;
% 
% for i=1:length(singularities)
%     % not using the edge indexing of shortest path, as I'm concerned about
%     % index issues
%     primalShortestPath = shortestpath(NBprimalGraph,singularities(i),closestVtx(i));
%     pathLength = length(primalShortestPath)-1;
%     
%     % convert to primal indices into top-left corner grid
%     [~, colPSP] = ind2sub([M+1, N+1], primalShortestPath);
%     TLprimalShortestPath = primalShortestPath - (colPSP - 1);
%     
%     % get edges in terms of primal vertex index endpoints (TL)
%     TLprimalEdges = [TLprimalShortestPath', circshift(TLprimalShortestPath,-1)'];
%     TLprimalEdges = TLprimalEdges(1:pathLength,:);
%     
%     % find locations in intEdges and adjust d_theta according to whether it
%     % agrees or disagrees with standard orientation
%     [~,addLocs] = ismember(TLprimalEdges,intEdgesOrig,'rows');
%     if max(addLocs) > 0
%         addLocs = nonzeros(addLocs);
%         d_theta(addLocs) = d_theta(addLocs) - singIndex(i)*2*pi;
%     end
%     [~,subLocs] = ismember(circshift(TLprimalEdges,1,2),intEdgesOrig,'rows');
%     if max(subLocs) > 0
%         subLocs = nonzeros(subLocs);
%         d_theta(subLocs) = d_theta(subLocs) + singIndex(i)*2*pi;
%     end
% end

%% Comb according to assignment and correct boundary indices

%ext 
OldIndices(1: numExteriorBoundaries) = 1 - H(1: (numExteriorBoundaries), :)*d_theta ./ (2 * pi);
%int 
OldIndices((numExteriorBoundaries + 1): numBoundaries) = -1 - H((numExteriorBoundaries + 1): numBoundaries, :)*d_theta ./ (2 * pi);

preclean_dtheta = d_theta;

singVisited = zeros(length(singularities),1);
plusCtr = 0;
minusCtr = 0;

% finalA=[[d0sub]';H] ;
% 
% b_1=zeros(size(d0sub,2),1);
% 
% finalb=[b_1;b];

%% Visualize singularities & indices

% scatter(singsxy(:,2)-0.5,singsxy(:,1)-0.5,40,'red')
% text(singsxy(:,2)-0.5,singsxy(:,1)-0.5,num2str(singIndex),'VerticalAlignment','bottom','HorizontalAlignment','right')
% set(gca,'YDir','reverse');
% hold on 

%%
for i=1:length(singularities)
    % skip if we've been there, but don't forget to update the counter
    if singVisited(i) == 1
        if singIndex(i) < 0
            minusCtr = minusCtr+1;
        else
            plusCtr = plusCtr+1;
        end
        continue
    end
    % check if it's a plus or minus
    if singIndex(i) < 0
        minusCtr = minusCtr+1;
        destination = find(assn(minusCtr,:));
        % see if it's going to a boundary or not
        if destination > length(singsPlus)
            % First, update index
            OldIndices(IMinus(minusCtr)) = OldIndices(IMinus(minusCtr)) + singIndex(i);

            % Now, let's adjust dtheta
            % find closest boundary vertex
            closestBdyVtx = closestMinusVtx(minusCtr,IMinus(minusCtr));
            
            % not using the edge indexing of shortest path, as I'm concerned about index issues
            primalShortestPath = shortestpath(NBminusPrimalDigraph,singularities(i),closestBdyVtx);
            pathLength = length(primalShortestPath)-1;
            
            % convert to primal indices into top-left corner grid
            [~, colPSP] = ind2sub([M+1, N+1], primalShortestPath);
            TLprimalShortestPath = primalShortestPath - (colPSP - 1);
            
            % get edges in terms of primal vertex index endpoints (TL)
            TLprimalEdges = [TLprimalShortestPath', circshift(TLprimalShortestPath,-1)'];
            TLprimalEdges = TLprimalEdges(1:pathLength,:);

            % find locations in intEdges and adjust d_theta according to whether it
            % agrees or disagrees with standard orientation
            [~,addLocs] = ismember(TLprimalEdges,intEdgesOrig,'rows');
            if max(addLocs) > 0
                addLocs = nonzeros(addLocs);
                d_theta(addLocs) = d_theta(addLocs) - singIndex(i)*2*pi;
            end
            [~,subLocs] = ismember(circshift(TLprimalEdges,1,2),intEdgesOrig,'rows');
            if max(subLocs) > 0
                subLocs = nonzeros(subLocs);
                d_theta(subLocs) = d_theta(subLocs) + singIndex(i)*2*pi;
            end
            
            % Check that boundary indices are satisfied
            % Form rotational constraints for boundary
            if IMinus(minusCtr) <= numExteriorBoundaries
                lineIndex = (OldIndices(IMinus(minusCtr)) + 1)/2;
                rotConstraint = (1-lineIndex)*2*pi;
            else
                lineIndex = (OldIndices(IMinus(minusCtr)) - 1)/2;
                rotConstraint = (-1-lineIndex)*2*pi;
            end
            if abs(H(IMinus(minusCtr),:)*(d_theta/2) - rotConstraint) > 1e-12
                assert(1 > 0);
            end
%% visualization            
%             [TLprimalEdgesrow, TLprimalEdgescol]= ind2sub(size(narrowBand),TLprimalEdges);
%             plot(TLprimalEdgescol-0.5, TLprimalEdgesrow-0.5, 'LineWidth', 3);
%             hold on
        else % if it maps to another vertex, let's find it 
            % find destination vertex & note we've been there
            destinationVtx = singularities(posSingsToSings(destination));
            singVisited(posSingsToSings(destination),1) = 1;
            
            % not using the edge indexing of shortest path, as I'm concerned about index issues
            primalShortestPath = shortestpath(NBminusPrimalDigraph,singularities(i),destinationVtx);
            pathLength = length(primalShortestPath)-1;
            
            % convert to primal indices into top-left corner grid
            [~, colPSP] = ind2sub([M+1, N+1], primalShortestPath);
            TLprimalShortestPath = primalShortestPath - (colPSP - 1);
            
            % get edges in terms of primal vertex index endpoints (TL)
            TLprimalEdges = [TLprimalShortestPath', circshift(TLprimalShortestPath,-1)'];
            TLprimalEdges = TLprimalEdges(1:pathLength,:);

            % find locations in intEdges and adjust d_theta according to whether it
            % agrees or disagrees with standard orientation
            [~,addLocs] = ismember(TLprimalEdges,intEdgesOrig,'rows');
            if max(addLocs) > 0
                addLocs = nonzeros(addLocs);
                d_theta(addLocs) = d_theta(addLocs) - singIndex(i)*2*pi;
            end
            [~,subLocs] = ismember(circshift(TLprimalEdges,1,2),intEdgesOrig,'rows');
            if max(subLocs) > 0
                subLocs = nonzeros(subLocs);
                d_theta(subLocs) = d_theta(subLocs) + singIndex(i)*2*pi;
            end
            %% visualization
%             [TLprimalEdgesrow, TLprimalEdgescol]= ind2sub(size(narrowBand),TLprimalEdges);
%             plot(TLprimalEdgescol-0.5, TLprimalEdgesrow-0.5, 'LineWidth', 3);
%             hold on
        end
                
    else % A plus singularity
        plusCtr = plusCtr+1;
        destination = find(assn(:,plusCtr));
        % see if it's going to a boundary or not
        if destination > length(singsMinus)
            % First, update index
            OldIndices(IPlus(plusCtr)) = OldIndices(IPlus(plusCtr)) + singIndex(i);

            % Now, let's adjust dtheta
            % find closest boundary vertex
            closestBdyVtx = closestPlusVtx(plusCtr,IPlus(plusCtr));
            
            % not using the edge indexing of shortest path, as I'm concerned about index issues
            primalShortestPath = shortestpath(NBplusPrimalDigraph,singularities(i),closestBdyVtx);
            pathLength = length(primalShortestPath)-1;
            
            % convert to primal indices into top-left corner grid
            [~, colPSP] = ind2sub([M+1, N+1], primalShortestPath);
            TLprimalShortestPath = primalShortestPath - (colPSP - 1);
            
            % get edges in terms of primal vertex index endpoints (TL)
            TLprimalEdges = [TLprimalShortestPath', circshift(TLprimalShortestPath,-1)'];
            TLprimalEdges = TLprimalEdges(1:pathLength,:);

            % find locations in intEdges and adjust d_theta according to whether it
            % agrees or disagrees with standard orientation
            [~,addLocs] = ismember(TLprimalEdges,intEdgesOrig,'rows');
            if max(addLocs) > 0
                addLocs = nonzeros(addLocs);
                d_theta(addLocs) = d_theta(addLocs) - singIndex(i)*2*pi;
            end
            [~,subLocs] = ismember(circshift(TLprimalEdges,1,2),intEdgesOrig,'rows');
            if max(subLocs) > 0
                subLocs = nonzeros(subLocs);
                d_theta(subLocs) = d_theta(subLocs) + singIndex(i)*2*pi;
            end
            
            % Check that boundary indices are satisfied
            % Form rotational constraints for boundary
            if IPlus(plusCtr) <= numExteriorBoundaries
                lineIndex = (OldIndices(IPlus(plusCtr)) + 1)/2;
                rotConstraint = (1-lineIndex)*2*pi;
            else
                lineIndex = (OldIndices(IPlus(plusCtr)) - 1)/2;
                rotConstraint = (-1-lineIndex)*2*pi;
            end
            if abs(H(IPlus(plusCtr),:)*(d_theta/2) - rotConstraint) > 1e-12
                assert(1 > 0);
            end
            %% visualization 
%             [TLprimalEdgesrow, TLprimalEdgescol]= ind2sub(size(narrowBand),TLprimalEdges);
%             plot(TLprimalEdgescol-0.5, TLprimalEdgesrow-0.5, 'LineWidth', 3);
%             hold on
        else % if it maps to another vertex, let's find it 
            % find destination vertex & note we've been there
            destinationVtx = singularities(negSingsToSings(destination));
            singVisited(negSingsToSings(destination)) = 1;
            
            % not using the edge indexing of shortest path, as I'm concerned about index issues
            primalShortestPath = shortestpath(NBplusPrimalDigraph,singularities(i),destinationVtx);
            pathLength = length(primalShortestPath)-1;
            
            % convert to primal indices into top-left corner grid
            [~, colPSP] = ind2sub([M+1, N+1], primalShortestPath);
            TLprimalShortestPath = primalShortestPath - (colPSP - 1);
            
            % get edges in terms of primal vertex index endpoints (TL)
            TLprimalEdges = [TLprimalShortestPath', circshift(TLprimalShortestPath,-1)'];
            TLprimalEdges = TLprimalEdges(1:pathLength,:);

            % find locations in intEdges and adjust d_theta according to whether it
            % agrees or disagrees with standard orientation
            [~,addLocs] = ismember(TLprimalEdges,intEdgesOrig,'rows');
            if max(addLocs) > 0
                addLocs = nonzeros(addLocs);
                d_theta(addLocs) = d_theta(addLocs) - singIndex(i)*2*pi;
            end
            [~,subLocs] = ismember(circshift(TLprimalEdges,1,2),intEdgesOrig,'rows');
            if max(subLocs) > 0
                subLocs = nonzeros(subLocs);
                d_theta(subLocs) = d_theta(subLocs) + singIndex(i)*2*pi;
            end
            [TLprimalEdgesrow, TLprimalEdgescol]= ind2sub(size(narrowBand),TLprimalEdges);
            %% visualization
%             plot(TLprimalEdgescol-0.5, TLprimalEdgesrow-0.5, 'LineWidth', 3);
%             hold on
        end
    end
    % Check that our updated conditions hold?
    
end

end
