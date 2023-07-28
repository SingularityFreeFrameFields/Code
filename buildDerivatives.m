function [d0sub, d1sub, narrowBandEdges, intEdgeInd, d1,intEdges,intVertices,narrowBandEdgesFull] = buildDerivatives(narrowBand)
%% Docstring
%    [d0sub, d1sub, narrowBandEdges, intEdgeInd] = buildDerivative(narrowBand) 
%    returns two sparse matricies d0sub and d1sub and two full matricies.
%    
%    narrowBandEdges are primal EDGES of the narrowBand with without any 
%    duplicates in [source, target] form with size Q x 2 (where Q is the 
%    the number of unique primal edges.
%
%    intEdgeInd are the row indicies that will find the INTERIOR edges in the
%    narrowBand.
%
%    The matrix d0sub is a |interiorEdges| x |interiorVerticies| 
%    matrix with the following rule: 
%        d0(i,j) = 1 if vertex j the SOURCE of edge i.
%        d0(i,j) == -1 if vertiex j is the TARGET of edge i.
%        d0(i,j) == 0 otherwise
%     
%    The matrix d1 is a |interior Faces| x |interiorEdges| with the 
%    following rule:
%        d1(i,j) = 1 if edge j is oriented COUNTER-CLOCKWISE on face i 
%        d1(i,j) = -1 if edge j is oriented CLOCKWISE on face i. 
%        d1(i,j) == 0 otherwise  
%    INPUT:
%        narrowBand: an M x N logical array
%
%    OUPUT:
%        d0sub: A sparse |intEdges| x |intVerticies| discrete GRADIENT operator
%        d1sub: A sparse |Faces| x |intEdges|  discrete CURL operator
%  
%   NOTE: |Verticies| and |Edges| are the number of primal verticies and
%   primal edges surrounding each pixel in the narrowBand.
%
%% A WORD ON INDEXING CONVENTIONS IN THIS CODE
%
% The narrowBand is an M x N logical array with 0's denoting white pixels
% and 1's denoting black pixels.
%
% The DUAL VERTICIES (pixels), are linearly ordered according to MATLAB's
% array indexing. The number of dual verticies is precisely the same as
% the number of faces (pixels).
%
% Consider the PRIMAL MESH of the narrowBand (the matrix that assigns a
% vertex to each corner of a pixel). This verticies of this mesh have a
% linear order according to the MATLAB array indexing.
%
% The VERTICIES array is |V| x 1 array of the LINEAR INDICIES of the
% primal verticies on the corners of each black pixel.
%
% The EDGES array is a |E| x 2 array where each row is of the form 
% [source, target] where both source and target are LINEAR INDICIES from
% VERTICIES.
%
%
% narrowBandVerticies and narrowBandEdges have an implicit ordering. For
% example, the 'n-th' vertex in the narrowBand is:
%
%     narrowBandVerticies(n) = vertex  
%
%  and the 'k-th' edge in the narrowBand is:
%
%     narrowBandEdges(k,:) = [source, target]
%
% intVertexInd are the indicies that will find the INTERIOR verticies in
% the narrowBand. For simplicity, narrowBandVerticies(intVertexInd) 
% returns a |interiorVerticies| x 1 array. The column vector is of the form
% [v1, v2 , v3, ..., vI, ...]' where vI is the primal vertex index of the
% i-th interior vertex.
%
% intEdgeInd are the row indicies that will find the INTERIOR edges in the
% narrowBand. For example, narrowBandEdges(intEdgeInd,:) is a 
% |interiorEdges| x 2 array where each row is of the form 
% [v_source, v_target] and v_source and v_target are primal indicies as in 
% intVertexInd.
%
%% Create Dual Vertex List

% dualV gives the linear indicies of the black pixels
dualV = find(narrowBand == 1);
numFaces = length(dualV);      
[M,~] = size(narrowBand);

% pad narrowBand and get corresponding pixel 
narrowBandFull = [narrowBand; zeros(1,size(narrowBand,2))];
dualVFull = find(narrowBandFull == 1);

%% Build Edge and Vertex Matricies

% FUNDAMENTAL FLAW BELOW: THE INCREMENTING AMOUNTS ASSUME AN M-BY-N GRID,
% WHILE THE RANGE ASSUMES AN (M+1)-BY-(N+1) GRID. IN PARTICULAR, FOR PIXELS
% ON THE BOTTOM OF THE IMAGE, THE BOTTOM THREE EDGES DON'T MAKE SENSE

% Allocate storage for the edges of the narrow band
% NOTE: |edges| <= 4 * |faces| and |verticies| <= 4*|faces|
Verticies = zeros(4 * numFaces, 1);
Edges = zeros(4*numFaces, 2);
EdgesFull = zeros(4*numFaces, 2);

for i = 1:numFaces
    j = dualV(i);    %i-th face (pixel) in the narrow band
    k = dualVFull(i);
    % edges are in the form [source, target]
    faceEdges = [j, j + M;                          % top edge
                 j+1, j;                            % left edge
                 j+1, j + (M + 1);                  % bottom edge
                 j + (M + 1), j + M];               % right edge
    % add 4 edges to the next 4 entries
    Edges(4 * (i-1) + 1: 4*(i-1) + 4, :) = faceEdges;
    
    % for "full" grid
    faceEdgesFull = [k, k + (M + 1);                    % top edge
                 k+1, k;                            % left edge
                 k+1, k + (M + 2);                  % bottom edge
                 k + (M + 2), k + (M + 1)];         % right edge
    % add 4 edges to the next 4 entries
    EdgesFull(4 * (i-1) + 1: 4*(i-1) + 4, :) = faceEdgesFull;
    % add 4 verticies to next 4 entries
    Verticies(4 * (i-1) + 1: 4*(i-1) + 4, :) = [j, j+1, j + M, j + M + 1];
end

% For primal mesh construction only (currently)
narrowBandEdgesFull = unique(EdgesFull, 'rows', 'first');

%% Find Interior Edges Indicies

% eliminate duplicate edges and the number of occurences of each row
[narrowBandEdges, ~, edgeOccurences] = unique(Edges, 'rows', 'first');
numEdges = length(narrowBandEdges);

% find the edges that are repeated exactly 2 times
repeatedEdgeInd = grouptransform(edgeOccurences, edgeOccurences, @numel) == 2;

% Of the edges that are repeated twice, remove the copies 
intEdges = unique(Edges(repeatedEdgeInd, :), 'rows');

% Linear Indicies of the interior edges of the narrowband with respect to
% the total edges of the narrowband
intEdgeInd = find(ismember(narrowBandEdges, intEdges, 'rows'));

%% Find Interior Verticies Indicies

% eliminate duplicates in Verticies list and obtain a list of duplicates
[narrowBandVerticies, ~, vertexOccurences] = unique(Verticies);
numVerticies = length(narrowBandVerticies);

% find the vertex indicies that are repeated exactly 4 times 
repeatedVertexInd = groupcounts(narrowBandVerticies(vertexOccurences)) == 4;

% obtain the verticies in the narrow band that are repeated 4 times
intVertices = narrowBandVerticies(repeatedVertexInd);

% Linear Indicies of the interior verticies of the narrowband with respect to
% the total verticies of the narrowband i.e narrowBandverticies(intEdgeInd)
% == list of interior verticies 
intVertexInd = find(ismember(narrowBandVerticies, intVertices));


%% Build d0 Matrix

d0 = spalloc(numEdges, numVerticies, numEdges * 2);

% Backwards array orders the primal verticies of the narrow band
backwards(narrowBandVerticies) = 1:numVerticies;

% Narrow band ordering of the source and targets
SourceInd = backwards(narrowBandEdges(:,1));
TargetInd = backwards(narrowBandEdges(:,2));

indTarg = [];
indSource = [];

for i = 1:numEdges
    indTarg = [indTarg; sub2ind(size(d0), i, TargetInd(i))];
    indSource = [indSource; sub2ind(size(d0), i, SourceInd(i))];
end

d0(indTarg) = -1;
d0(indSource) = 1;

%% Build d1 Matrix

d1 = sparse(numFaces, numEdges, numFaces * 4);
% largest primal vertex index on narrow band (according to global
% coordinates

maxV = max(narrowBandVerticies);

% For [source target], indEdge returns the linear index of (source, target)
% which can be used in backwardEdges

indEdges = sub2ind([maxV maxV], narrowBandEdges(:,1), narrowBandEdges(:,2));

% backwardEdges orders the edges of the narrow band, i.e 
% backwardEdges(i,j) = k if (i j) is the k-th edge in the narrow band

backwardsEdges = sparse(maxV, maxV);
backwardsEdges(indEdges) = 1:numEdges;   % apparently this line takes a long time

faceEdges =  [dualV, dualV+1, dualV+1, dualV + (M + 1), dualV + M,  dualV,  dualV + (M + 1),  dualV + M];

posEdges = sub2ind([maxV maxV], faceEdges(:,3:4), faceEdges(:,7:8))';
negEdges = sub2ind([maxV maxV], faceEdges(:,1:2), faceEdges(:,5:6))';
PosFaceEdgeEntriesI = [1:numFaces; 1:numFaces];
PosFaceEdgeEntriesJ = backwardsEdges(posEdges);
NegFaceEdgeEntriesI = PosFaceEdgeEntriesI;
NegFaceEdgeEntriesJ = backwardsEdges(negEdges);

nPos = numel(PosFaceEdgeEntriesI);
nNeg = numel(NegFaceEdgeEntriesI);
d1 = sparse([PosFaceEdgeEntriesI(:); NegFaceEdgeEntriesI(:)], [PosFaceEdgeEntriesJ(:); NegFaceEdgeEntriesJ(:)], [ones(nPos,1); -ones(nNeg,1)]);
% Obtain only the sub-matricies of d0 and d1 that correspond to interior
% Vertex and Edges
d0sub = d0(intEdgeInd, intVertexInd);
d1sub = d1(:, intEdgeInd);
