function [primalGraph, NBprimalGraph] = PrimalMesh(narrowBand, narrowBandEdges, narrowBandEdgesFull)
%% Docstring
%    [narrowBandprimalGraph] = PrimalMesh(narrowBand) returns a 
%    the primal graph of a narrowBand (B/W image).
%
%    The Nodes of narrowBandprimalGraph correspond to the corners of 
%    black pixels (pixels that are turned on).
%    The Edges of narrowBandprimalGraph correspond to the sides of all
%    black pixels.
%
%    INPUT:
%        narrowBand: an M x N logical array
%
%    OUPUT:
%        narrowBandprimalGraph: a 'graph' object
% 

%% Node Construction

[M,N] = size(narrowBand);

% index each corner of each pixel in column - major order
% NOTICE: we use M+1 because we must account for the pixels at the
% bottom of the image. Furthermore, we use only use N because we want to
% purposefully omit the last column on the right in the next line of 
% code.
i = 1:N*(M+1);

% each pixel with index i, is connected via horizontal edge to pixel 
% i + (M+1) (the pixel in the column to the right)
j = i + (M+1)';

%% Edge Construction

% Orientation of Primal Edges:
%
%    ^>>>>>^
%    ^ ...............^
%    ^ ...............^
%    ^>>>>>^
%

% Horizontal Edges
H_edges = [i',j'];

% Vertical Edges

% index all (N+1) * (M + 1) pixel corners this time
k = 1:(N+1)*(M+1);
top_verts = 1:(M+1):((N)*(M+1) + 1);
% remove the pixel corners on the top row
s = setdiff(k, top_verts);
% each pixel in s is connected via vertical edge to pixel s - 1 (vertex above)
p = s - 1;
V_edges = [s',p'];
Edges = [H_edges; V_edges];

%% "Background" Primal Graph Construction

primalGraph = graph(Edges(:,1), Edges(:,2));
primalGraph.Nodes.PrimalIndex = (1: (N+1)*(M+1))';
primalGraph.Edges.PrimalIndex = (1:length(Edges))';

%% Construct narrowBand (primal) subgraph

% [~, C] = ind2sub([M,N], narrowBandEdges);
% tNarrowBandEdges = narrowBandEdges + (C-1);

NBprimalGraph = graph;
NBprimalGraph = addnode(NBprimalGraph, primalGraph.Nodes);
NBprimalGraph = addedge(NBprimalGraph, narrowBandEdgesFull(:,1), narrowBandEdgesFull(:,2));

[~, C] = ind2sub([M,N], narrowBandEdges);
tNarrowBandEdges = narrowBandEdges + (C-1);

NBprimalGraph2 = graph;
NBprimalGraph2 = addnode(NBprimalGraph2, primalGraph.Nodes);
NBprimalGraph2 = addedge(NBprimalGraph2, tNarrowBandEdges(:,1), tNarrowBandEdges(:,2));

% Uncomment to plot the graph
%{ 
 figure;
 hold on;
 h = plot(primalGraph);
 highlight(h, NBprimalGraph.Edges.EndNodes(:,1), NBprimalGraph.Edges.EndNodes(:,2), 'EdgeColor', 'black', 'NodeColor', 'red', 'LineWidth', 5);
%}
end

