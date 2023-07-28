tic

for filename=["my_inputs/narrowJunction.png"] %INPUT: png/jpeg images
        
close all; clearvars -except filename;

%config
runMatlabOptimization = false; %select 'true' to run the nonconvex optimization using Matlab's fmincon
runCppOptimization = true; %select 'true' to run the nonconvex optimization using IPOPT (C++)
loadCppOptimization = true; %must be 'true' if the C++ optimization is run 

%visualizations 
visualizeTangentAssignment = false;
visualizeInitialField = false;
visualizeBoundaries = false;
visualizeGradient = false;

%parameters
alignment_wt_opt = 2; %alignment term weight for the nonconvex optimization 
barrier_wt_opt = 1; %barrier term weight for the nonconvex optimization
smoothness_wt_opt = 10; %smoothness term weight for the nonconvex optimization
alignment_wt_init = 1; %alignment term weight for the convex optimization 
smoothness_wt_init = 8; %smoothness term weight for the nonconvex optimization 
wDiff = 0.0;
narrowBandIntensityThreshold = 90;
antialign = 0.05;

%Run the Vectorization algorithm of [Bessmeltsev and Solomon 2019] with singularity-free frame fields
pathToPolyVector = 'sgi_vectorization_polyvector/polyvector_thing.exe';

img = imread(filename);

bwImg = 255-(img(:,:,1)*0.299+img(:,:,2)*0.587+img(:,:,3)*0.114);%convert image to greyscale
m = size(bwImg,1);
n = size(bwImg,2);
%pad image on bottom and right to make sure all interior edges may be expressed
bwImg = [bwImg; zeros(1,n)]; %pad below
bwImg = [bwImg, zeros(m+1,1)]; %pad right
m = m + 1; %adjust m
n = n + 1; %adjust n

narrowBand_init = bwImg > narrowBandIntensityThreshold;
narrowBand = narrowBand_init;
for k=1:3
    narrowBand = repairMask(narrowBand);
end

[M,~] = size(narrowBand);
[Gmag, Gdir] = imgradient(bwImg);

%%
initializationWorked = false;

while ~initializationWorked
    [BoundaryPixelInd, numBoundaries, N, smallLenBoundaryInd, BoundaryPixelIndold, Nold, narrowBand]=indexBoundaryPixelsEd(M, narrowBand);
    [d0sub, d1sub, narrowBandEdges, intEdgeInd, d1,intEdges,intVerticies, narrowBandEdgesFull] = buildDerivatives(narrowBand);
    [H, dualedgeInd]=finaledge(M, BoundaryPixelInd, intEdgeInd, narrowBandEdges, numBoundaries);
    
    %needs Image Graphs toolbox: https://www.mathworks.com/matlabcentral/fileexchange/53614-image-graphs
    G = binaryImageGraph(narrowBand,4);

    nDualEdges=size(dualedgeInd,1);
    nPixelsNarrowBand = size(G.Nodes,1);

    [Tree, pred] = minspantree(G,'Type', 'forest');
    edges_bfs = bfsearch(Tree,1,'edgetonew','Restart',true);

    pixelind_bfs = bfsearch(Tree,1,'Restart',true); %pixelind_bfs is the list of graph vertex ids in the BFS order

    tree_edges = table2array(Tree.Edges);
    roots_bfs = [];
    connectedComponents = conncomp(Tree);
    nComponents = max(connectedComponents);

    myBFSorder = zeros(1,nPixelsNarrowBand);%for a vertex ID, myBFSorder returns its BFS order, i.e. it's the inverse of pixelind_bfs
    myBFSorder(pixelind_bfs) = 1:nPixelsNarrowBand;
    for i=1:nComponents
        myConnectedComponent = find(connectedComponents == i);
        [~,argmin] = min(myBFSorder(myConnectedComponent));
        roots_bfs(i) = myConnectedComponent(argmin);
    end

    %normalization factors 
    alignment_factor = nnz(Gmag);
    barrier_factor = nPixelsNarrowBand;
    smoothness_factor =nDualEdges;

    [primalGraph, NBprimalGraph] = PrimalMesh(narrowBand, narrowBandEdges,narrowBandEdgesFull);
    Initialization = FieldInitialization(d1, d1sub, bwImg, narrowBand, alignment_factor, smoothness_factor,...
        alignment_wt_init, smoothness_wt_init);
   
    if visualizeInitialField
    figure
    VisualizeInitialization(Initialization, narrowBand);
    figure
    set(gca,'YDir','reverse');
    hold on
    VisualizeInitialization(Initialization, narrowBand);
    end
    
    [Indices,initial_dtheta,singsxy,singIndex,~,backwardsLinearNarrowBand,preclean_dtheta] = CleanInitialization(Initialization, narrowBand,...
            d0sub, H, intEdges, primalGraph, NBprimalGraph, intVerticies, numBoundaries, N, smallLenBoundaryInd, BoundaryPixelInd, Nold);
    initializationWorked = true;
end

%convert cross field index to line field index (Eq 9) 
index_1=[(Indices(1:N) + 1 )/2]';
index_2=[(Indices((N+1):numBoundaries) - 1)/2]';

disp(['Exterior index: ',num2str(index_1')])
disp(['Interior indices: ', num2str(index_2')])

b=righthandside_b(numBoundaries, N, index_1, index_2);
[solangles, finalA, finalb]=solvingAngle(H, d0sub, b, initial_dtheta/2);

tau = exp(1i*(Gdir*pi/180+pi/2));
g = exp(1i*(Gdir*pi/180));
tau(Gmag < max(max(Gmag))./10) = 0;
Gmag(Gmag < max(max(Gmag))/10) = 0;
Gdir = deg2rad(Gdir);
tau2 = tau.^2;
tau2re = real(tau2);
tau2im = imag(tau2);

[a,b] = ndgrid(1:m,1:n);
u=Gmag.*cos(Gdir);
v=Gmag.*sin(Gdir);
u(~bwImg)=0;
v(~bwImg)=0;

%% visualize gradient
if visualizeGradient
q1 = quiver(b(:),a(:),u(:),-v(:),'r');
set(gca,'YDir','reverse');
hold on;
q1.ShowArrowHead = 'on';
axis equal;
end

%% visualize boundaries
if visualizeBoundaries
figure
set(gca,'YDir','reverse');
hold on;
axis equal;
for n=1:numel(BoundaryPixelInd)
    [i,j] = ind2sub(size(narrowBand),BoundaryPixelInd{n});
    lineObject = plot(j,i);
    color = get(lineObject, 'Color');
    centerIndex = round(numel(i/2));
    centralPt = [i(centerIndex); j(centerIndex)];
    text(centralPt(2), centralPt(1),strcat("#",int2str(n)),'Color',color);
end
end

%% optimization set-up
[ddtheta2_dtheta2,alledges_bfs] = theta2_derivatives(pixelind_bfs, G,dualedgeInd, narrowBand);
[dtheta1_ddtheta1,ddtheta2_ddtheta1] =  theta1_derivatives(roots_bfs, edges_bfs, pixelind_bfs, G, dualedgeInd);
g_theta_BFS = 2*Gdir(Tree.Nodes.PixelIndex(pixelind_bfs));
g_magni_BFS = Gmag(Tree.Nodes.PixelIndex(pixelind_bfs));
ind_of_zero_gradients_BFS = find(g_magni_BFS==0);

%% normalized optimization weights
wA = alignment_wt_opt/alignment_factor;
wB = barrier_wt_opt/barrier_factor;
wS = smoothness_wt_opt/smoothness_factor;

data.edges_bfs = edges_bfs;
data.roots_bfs = roots_bfs;
data.pixelind_bfs = pixelind_bfs;
data.ind_of_zeros = ind_of_zero_gradients_BFS;
data.Tree = Tree;
data.dualedgeInd = dualedgeInd;
data.alledges_bfs = alledges_bfs;
data.narrowBand = narrowBand;
data.myBFSorder = myBFSorder;
data.G = G;
data.g_theta = g_theta_BFS';
data.nDualEdges = nDualEdges;
data.nComponents = nComponents;
data.antialign = antialign;
phi_init = angle(Initialization(roots_bfs)+1i*Initialization(roots_bfs+nPixelsNarrowBand));
init = [phi_init/2; initial_dtheta/2; ones(nPixelsNarrowBand,1)*pi];
data.init = init;

%% matlab optimization set-up
fun = @(x) myfunction_v( x(1:nComponents), x((nComponents+1):(nDualEdges+nComponents)), x(nDualEdges+1+nComponents:end),...
    wA, wB, wS, wDiff, ddtheta2_dtheta2, dtheta1_ddtheta1,ddtheta2_ddtheta1, data);
funHessian = @(x,lambdas) myfunction_hessian( x(1:nComponents), x((nComponents+1):(nDualEdges+nComponents)), x(nDualEdges+1+nComponents:end),lambdas,...
    wA, wB, wS, wDiff, ddtheta2_dtheta2, dtheta1_ddtheta1,ddtheta2_ddtheta1, data);
linConstraintRows = size(finalA,1);
finalA = [sparse(linConstraintRows,size(roots_bfs,2)) finalA sparse(linConstraintRows,nPixelsNarrowBand)];

%% export for C++
if runCppOptimization
exportFilename = strcat(filename,datestr(now, 'mmm_dd_HH-MM-SS'),'data.m');
export_for_cpp(exportFilename,phi_init/2, initial_dtheta/2, ones(nPixelsNarrowBand,1)*pi,wA, wB, wS, antialign, ddtheta2_dtheta2, dtheta1_ddtheta1, ddtheta2_ddtheta1, edges_bfs, roots_bfs, pixelind_bfs, ind_of_zero_gradients_BFS, Tree, dualedgeInd, alledges_bfs, myBFSorder, G, G.Nodes.PixelIndex, g_theta_BFS, finalA, finalb);
end

%% matlab optimization set-up (cont.)
if runMatlabOptimization
options = optimoptions('fmincon','Display','iter','MaxIterations',10000,'OptimalityTolerance',1e-6,'MaxFunctionEvaluations',3.000000e+10,'SpecifyObjectiveGradient',true,'HessianFcn',funHessian);
[d_theta1_and_theta2,fval,exitflag,output] = fmincon(fun, init, [], [], finalA, finalb, [], [], [], options);
d_theta1 = d_theta1_and_theta2((nComponents+1):(nDualEdges+nComponents));
phi = d_theta1_and_theta2(1:nComponents);
theta2 = d_theta1_and_theta2((nDualEdges+nComponents+1):end);
end

%% C++ optimization 
if runCppOptimization || loadCppOptimization
tic
if runCppOptimization
system( sprintf('\"%s\" %s -ff2 %s', pathToPolyVector, filename, exportFilename));
end
filename=convertStringsToChars(filename);
cppOutFilename = filename(1:end-4)+"_opt.m";
cppOutFilenameFixed = strrep(cppOutFilename,"-","_"); 
if runCppOptimization && (cppOutFilename~=cppOutFilenameFixed) 
    movefile(cppOutFilename,cppOutFilenameFixed);
end
run(cppOutFilenameFixed);
toc
end

%%
figure;
set(gca,'YDir','reverse');
hold on; axis equal;
[theta1_output] = fieldReconstruction(phi, d_theta1, theta2, narrowBand, edges_bfs, roots_bfs, pixelind_bfs, Tree, G, dualedgeInd); %output for vectorization testing

finalinfo_theta2 = table2array(Tree.Nodes);
angleforvec2 = zeros(size(narrowBand));
thirdcolumn2 = finalinfo_theta2(:,3);

for theta2_index=1:size(pixelind_bfs,1)
    angleforvec2(thirdcolumn2(pixelind_bfs(theta2_index)))=theta2(theta2_index);
end
theta2_output = [angleforvec2]; %output for vectorization testing
toc

%% tangent assignment
utau = Gmag.*cos(Gdir+pi/2);
vtau = Gmag.*sin(Gdir+pi/2);
dist1 = abs(exp(1i*theta1_output)-exp(2*1i*(Gdir+pi/2)));
dist2 = abs(exp(1i*(theta1_output+theta2_output))-exp(2*1i*(Gdir+pi/2))); 
dist1(~narrowBand) = inf;
dist2(~narrowBand) = inf;

if visualizeTangentAssignment
hold on;
quiver(b(:),a(:),utau(:),-vtau(:),'r');
figure
imagesc(dist1);
figure
imagesc(dist2); axis equal
end

Limg = laplacian(G);
gEdges = table2array(G.Edges);

graphD = sparse(nDualEdges, nPixelsNarrowBand);
graphD(sub2ind(size(graphD),(1:nDualEdges)',gEdges(:,1))) = -1;
graphD(sub2ind(size(graphD),(1:nDualEdges)',gEdges(:,2))) = 1;

dTheta1 = graphD*theta1_output(narrowBand);
jumpEdges = find(abs(dTheta1)>pi/2); %discontinuities of theta1 function over the graph edges -- do not enforce smoothness over them
graphD(jumpEdges,:) = []; %remove those

Lmine = graphD'*graphD; %weighted laplacian

Wimp = ones(nPixelsNarrowBand,1);
Wimp(Gmag(narrowBand)==0) = 0;
Wimp = spdiags(Wimp,0,nPixelsNarrowBand,nPixelsNarrowBand);
dist1Sm = lsqminnorm(5*Lmine+Wimp,Wimp*dist1(narrowBand));
dist2Sm = lsqminnorm(5*Lmine+Wimp,Wimp*dist2(narrowBand));
dist1Sm_full = ones(size(dist1))*inf;
dist2Sm_full = ones(size(dist2))*inf;
dist1Sm_full(narrowBand) = dist1Sm;
dist2Sm_full(narrowBand) = dist2Sm;

if visualizeTangentAssignment
figure
imagesc(dist1Sm_full)
figure
imagesc(dist2Sm_full)
end

tangentAssignment = exp(0.5*1i*(theta1_output+theta2_output));
theta1_mask = dist1Sm_full<dist2Sm_full;
tangentAssignment(theta1_mask) = exp(0.5*1i*theta1_output(theta1_mask));
tangentAssignment(~narrowBand) = 0;

if visualizeTangentAssignment
figure
imagesc(dist1Sm_full>dist2Sm_full)
axis equal
end

%%
exportFilename = strcat(filename,datestr(now, 'mmm_dd_HH-MM-SS'),'.m');
export_for_vectorization(exportFilename,narrowBand,theta1_output,theta2_output,tangentAssignment);
system( sprintf('\"%s\" %s -ff %s > vectorizationLog.txt', pathToPolyVector, filename, exportFilename));
end

