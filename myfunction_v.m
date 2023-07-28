function [fun, grad]= myfunction_v (phi, d_theta1, theta2, wA, wB, wS, wDiff, ddtheta2_dtheta2, dtheta1_ddtheta1,ddtheta2_ddtheta1, q)

%NOTE: theta1, theta2, and g_theta are all in BFS order
[theta1, ~] = integrate_theta(phi, d_theta1, q.edges_bfs, q.roots_bfs, q.pixelind_bfs, q.Tree, q.G, q.dualedgeInd, q.narrowBand);

for nr_dualedge=1:size(q.dualedgeInd,1)%-size(roots_bfs,2)
    d_theta2(nr_dualedge) = theta2(q.pixelind_bfs==q.alledges_bfs(nr_dualedge,2)) - theta2(q.pixelind_bfs==q.alledges_bfs(nr_dualedge,1)) + d_theta1(nr_dualedge);
end

%only for non-zero values
%g_theta = 2*Gdir(Tree.Nodes.PixelIndex(pixelind_bfs));
%g_magni = Gmag(Tree.Nodes.PixelIndex(pixelind_bfs));
%ind_of_zeros = find(g_magni==0);
theta1_small=theta1;
g_theta = q.g_theta';
g_thetasmall = g_theta;
theta2_small=theta2;
theta1_small(q.ind_of_zeros) = [];
theta2_small(q.ind_of_zeros) = [];
g_thetasmall(q.ind_of_zeros) = [];
antialign = q.antialign;

alignment = sum(antialign.*(-cos(theta1_small-g_thetasmall)+1).*(-cos(theta1_small+theta2_small-g_thetasmall)+1) + (cos(theta1_small-g_thetasmall)+1).*(cos(theta1_small+theta2_small-g_thetasmall)+1));

barrier = sum(-log(theta2)-log(2*pi-theta2) + 2*log(pi));
if any((theta2<0) | (theta2>2*pi))
    barrier = inf;
end

d_theta2=d_theta2';
smoothness = d_theta1'*d_theta1 + d_theta2'*d_theta2;

diff_init = [phi; d_theta1; theta2]-q.init;
fun = wA*alignment+ wB*barrier + wS*smoothness + 0.5*wDiff*(diff_init'*diff_init);

if nargout > 1
    nRoots = numel(q.roots_bfs);
    nPixels = numel(theta2);
    g_theta(q.ind_of_zeros)=theta1(q.ind_of_zeros);%pixel with gradient zero
    dtheta1_dphi = zeros(nPixels, nRoots);

    root_BFSorder = q.myBFSorder(q.roots_bfs);

    for i=1:nRoots-1
        dtheta1_dphi(root_BFSorder(i):root_BFSorder(i+1)-1,i)=1;
    end
    dtheta1_dphi(root_BFSorder(nRoots):end,nRoots)=1;
   
    dA_dtheta1_antialign = antialign.*(1-cos(theta1-g_theta)).*sin(theta1+theta2-g_theta)+antialign.*sin(theta1-g_theta).*(1-cos(theta1+theta2-g_theta));
    dA_dtheta1 = -((cos(theta1-g_theta)+1).*sin(theta1+theta2-g_theta))-(sin(theta1-g_theta).*(cos(theta1+theta2-g_theta)+1))+dA_dtheta1_antialign;
    dA_dtheta1(q.ind_of_zeros)=0; 

    dA_dphi = dtheta1_dphi'*dA_dtheta1;
    dA_ddtheta1 = dA_dtheta1'*dtheta1_ddtheta1;
    dA_dtheta2_antialign = antialign.*(1-cos(g_theta-theta1)).*sin(theta1+theta2-g_theta);
    dA_dtheta2_prelim = -((cos(g_theta-theta1)+1).*sin(theta1+theta2-g_theta))+dA_dtheta2_antialign;
    dA_dtheta2 = dA_dtheta2_prelim';
    dA_dtheta2(q.ind_of_zeros)=0;

    dbarrier_dphi = zeros(size(phi));
    dbarrier_ddtheta1 = zeros(size(d_theta1))';
    dbarrier_dtheta2 = [-1./theta2 + 1./(2*pi-theta2)]';

    dsmoothness_dphi = 0;
    dsmoothness_ddtheta1 = 2*d_theta1' + (2*d_theta2')*ddtheta2_ddtheta1;
    dsmoothness_dtheta2 = (2*d_theta2')*ddtheta2_dtheta2;

    dx_wrt_phi = wA*dA_dphi + wB*dbarrier_dphi + wS*dsmoothness_dphi;
    dx_wrt_dtheta1 =  wA*dA_ddtheta1 +wB*dbarrier_ddtheta1 + wS*dsmoothness_ddtheta1;
    dx_wrt_dtheta2 = wA*dA_dtheta2 +wB*dbarrier_dtheta2 + wS*dsmoothness_dtheta2;

    grad = [dx_wrt_phi' dx_wrt_dtheta1 dx_wrt_dtheta2]+wDiff*diff_init';
end
end
