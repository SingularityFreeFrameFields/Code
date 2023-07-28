function [BoundaryPixelInd, numBoundaries, NEight, smallLenBoundaryInd, BoundaryPixelIndold, Nold, narrowBand] = indexBoundaryPixels(M, narrowBand)
%% Docstring
%    [BoundaryPixelInd, numBoundaries, N] = indexBoundaryPixels(img) returns 
%    a N x 1 CELL where the i-th array in the cell is a list of pixels in 
%    the i-th boundary of the narrowBand.
%    Each list in the cell holds the LINEAR indicies of 
%    each pixel and are ordered such that EXTERIOR boundaries are
%    orientated COUNTER-CLOCKWISE and INTERIOR boundaries are orientated
%    CLOCKWISE.  
%   
%
%    INPUT: 
%        narrowBand: a BINARY array corresponding an image
%
%    OUTPUT:
%        BoundaryPixelInd: An N x 1 cell containing all of the linear indicies ofthe 
%        boundaries of narrowBand with the proper orientation for interior and 
%        exterior verticies.
%
%        numBoundaries: the total number of boundaries   

%        N: the number of exterior boundaries
close all;                      
%% Create Interior Boundary Edges
 
% 8-connectivity gives 4-connectivity for interior boundaries
[B,~,NEight,~] = bwboundaries(narrowBand, 8);

% Initialize a cell with an array for each boundary
BoundaryPixelInd = cell(size(B));

% The number of boundaries resulting from 8-connectivity and 4-connectivity
% are different.
BoundaryPixelIndInts = cell(size(B)-[NEight, 0]); %Just for interior boundaries

%        [B,L,N,A] = bwboundaries(narrowBand, 8);
%        figure; 
%        imshow(narrowBand); 
%        hold on;
%        colors = ['b' 'g' 'r' 'c' 'm' 'y'];
%        for k = 1:length(B)
%            boundary = B{k};
%            cidx = mod(k,length(colors))+1;
%            plot(boundary(:,2), boundary(:,1), colors(cidx), 'LineWidth',2);
%            % Randomize text position for better visibility
%            rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
%            col = boundary(rndRow,2); row = boundary(rndRow,1);
%            h = text(col+1, row-1, num2str(L(row,col)));
%            set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
%        end
%        figure;
%        spy(A);


for k = NEight+1:length(B)
    boundary = B{k};
    %BoundaryPixelInd{k} = sub2ind(size(narrowBand),boundary(:,1), boundary(:,2));
    BoundaryPixelIndInts{k-NEight} = sub2ind(size(narrowBand),boundary(:,1), boundary(:,2));
    %remove the boundary arrays that are very small and thus irrelevant
%     if size(BoundaryPixelInd{k},1) > 10
%         Diff = BoundaryPixelInd{k} - circshift(BoundaryPixelInd{k}, -1);
%         Consec = find(Diff== -M);
%         if length(Consec) == 0
%             startpt = BoundaryPixelInd{k}(1);
%         else
%             Horiedge = Consec(1);
%             startpt = BoundaryPixelInd{k}(Horiedge);
%         end
%         [r1,c1] = ind2sub(size(narrowBand),startpt);
%         % bwtraceboundary is used because bwboundaries incorrectly
%         % identifies interior boundaries
%         C = bwtraceboundary(narrowBand, [r1, c1-1], 'e', 4, inf, 'clockwise');
%         BoundaryPixelInd{k} = sub2ind(size(narrowBand),C(:,1), C(:,2));
%     end   

    % Processing interior boundaries separately
    %if size(BoundaryPixelIndInts{k-NEight},1) > 10
        Diff = BoundaryPixelIndInts{k-NEight} - circshift(BoundaryPixelIndInts{k-NEight}, -1);
        Consec = find(Diff== -M);
        if length(Consec) == 0
            startpt = BoundaryPixelIndInts{k-NEight}(1);
        else
            Horiedge = Consec(1);
            startpt = BoundaryPixelIndInts{k-NEight}(Horiedge);
        end
        [r1,c1] = ind2sub(size(narrowBand),startpt);
        % bwtraceboundary is used because bwboundaries incorrectly
        % identifies interior boundaries
        C = bwtraceboundary(narrowBand, [r1, c1-1], 'e', 4, inf, 'clockwise');
        BoundaryPixelIndInts{k-NEight} = sub2ind(size(narrowBand),C(:,1), C(:,2));
    %end   
end

%% Create Exterior Boundary Edges


[B,~,NFour,~] = bwboundaries(narrowBand,4);

BoundaryPixelIndExts = cell([NFour, 1]); %Just for exterior boundaries

% [B,L,N,~] = bwboundaries(narrowBand, 4);
%        figure; 
%        imshow(narrowBand); 
%        hold on;
%        colors = ['b' 'g' 'r' 'c' 'm' 'y'];
%        for k = 1:length(B)
%            boundary = B{k};
%            cidx = mod(k,length(colors))+1;
%            plot(boundary(:,2), boundary(:,1), colors(cidx), 'LineWidth',2);
%            % Randomize text position for better visibility
%            rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
%            col = boundary(rndRow,2); row = boundary(rndRow,1);
%            h = text(col+1, row-1, num2str(L(row,col)));
%            set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
%        end

for k = 1:NFour
    boundary = B{k};
    % make it ccw
    flip(boundary);
    BoundaryPixelIndExts{k} = sub2ind(size(narrowBand),boundary(:,1), boundary(:,2));
end  

% Gathers all 4-connected exterior boundaries, regardless of size; and all
% interior boundaries (in the narrowband and 4-connected if large enough,
% out of the narrowband and 4-connected if not)
BoundaryPixelInd = [BoundaryPixelIndExts;BoundaryPixelIndInts];

BoundaryPixelIndold =  BoundaryPixelInd;

%% boundary length threshold

[nrows,~]=cellfun(@size,BoundaryPixelInd,'UniformOutput',false);
smallLenBoundaryInd = find(cell2mat(nrows)<10);

% Eliminates small exterior boundaries
for i=1:length(smallLenBoundaryInd)
    narrowBand(BoundaryPixelInd{smallLenBoundaryInd(i)}) = 0;
end

BoundaryPixelInd(find(cell2mat(nrows)<10)) = [];
numBoundaries = length(BoundaryPixelInd);

%updating number of exterior boundaries
Nold = NEight;

goodLenBoundaryInd = find(~(cell2mat(nrows)<10));
NEight = numel(find(goodLenBoundaryInd<=NEight));


