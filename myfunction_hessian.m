function [H] = myfunction_hessian(phi, d_theta1, theta2, lambdas, wA, wB, wS, wDiff, ddtheta2_dtheta2, dtheta1_ddtheta1,ddtheta2_ddtheta1, q)

[theta1, ~] = integrate_theta(phi, d_theta1, q.edges_bfs, q.roots_bfs, q.pixelind_bfs, q.Tree, q.G, q.dualedgeInd, q.narrowBand);
nRoots = numel(q.roots_bfs);
nPixels = numel(theta2);
nDualEdges = numel(d_theta1);

nVars = numel(phi)+numel(d_theta1)+numel(theta2);

H = sparse(nVars, nVars);
g_theta = q.g_theta';
antialign = q.antialign;

for nr_dualedge=1:size(q.dualedgeInd,1)
    d_theta2(nr_dualedge) = theta2(q.pixelind_bfs==q.alledges_bfs(nr_dualedge,2)) - theta2(q.pixelind_bfs==q.alledges_bfs(nr_dualedge,1)) + d_theta1(nr_dualedge);
end

root_BFSorder = q.myBFSorder(q.roots_bfs);
dtheta1_dphi = zeros(nPixels, nRoots);
for i=1:nRoots-1
    dtheta1_dphi(root_BFSorder(i):root_BFSorder(i+1)-1,i)=1;
end
dtheta1_dphi(root_BFSorder(nRoots):end,nRoots)=1;

d2A_dtheta1_2_diag_antialign = 2*antialign*sin(theta1-g_theta).*sin(theta1-g_theta+theta2)+(antialign-2*antialign.*cos(theta1-g_theta)).*cos(theta1-g_theta+theta2)+antialign*cos(theta1-g_theta);
d2A_dtheta1_2_diag = - cos(theta1 - g_theta + theta2) - 2*cos(2*theta1 - 2*g_theta + theta2) - cos(g_theta - theta1) + d2A_dtheta1_2_diag_antialign;    
d2A_dtheta1_2_diag(q.ind_of_zeros) = 0;
d2A_dtheta1_2 = spdiags(d2A_dtheta1_2_diag,0,nPixels,nPixels);

d2A_dtheta1_dtheta2_diag_antialign =  antialign*(sin(theta1-g_theta).*sin(theta1-g_theta+theta2)+(1-cos(theta1-g_theta)).*cos(theta1-g_theta+theta2));
d2A_dtheta1_dtheta2_diag = - cos(theta1 - g_theta + theta2).*(cos(g_theta - theta1) + 1) - sin(g_theta - theta1).*sin(theta1 - g_theta + theta2) + d2A_dtheta1_dtheta2_diag_antialign;
d2A_dtheta1_dtheta2_diag(q.ind_of_zeros) = 0;
d2A_dtheta1_dtheta2=  spdiags(d2A_dtheta1_dtheta2_diag,0,nPixels,nPixels);

d2A_ddtheta1_2 = dtheta1_ddtheta1'*d2A_dtheta1_2*dtheta1_ddtheta1;%should be nDualEdges x nDualEdges
d2A_ddtheta1_dtheta2 = d2A_dtheta1_dtheta2*dtheta1_ddtheta1;
d2A_ddtheta1_dphi = dtheta1_dphi'*d2A_dtheta1_2*dtheta1_ddtheta1;
d2A_dtheta2_dphi = d2A_dtheta1_dtheta2*dtheta1_dphi;
d2A_dphi2 =  dtheta1_dphi'*d2A_dtheta1_2*dtheta1_dphi;

d2A_dtheta2_2_diag_antialign = antialign*(1-cos(g_theta-theta1)).*cos(theta2-g_theta+theta1);
d2A_dtheta2_2_diag = -cos(theta1 - g_theta + theta2).*(cos(g_theta - theta1) + 1) + d2A_dtheta2_2_diag_antialign;
d2A_dtheta2_2_diag(q.ind_of_zeros) = 0;
d2A_dtheta2_2 = spdiags(d2A_dtheta2_2_diag,0,nPixels,nPixels);

d2B_dtheta2_2_diag = 1./theta2.^2 + 1./(theta2 - 2*pi).^2;
d2B_dtheta2_2 = spdiags(d2B_dtheta2_2_diag,0,nPixels,nPixels);

d2S_dtheta1_2 = 2*speye(nDualEdges,nDualEdges)+2*(ddtheta2_ddtheta1'*ddtheta2_ddtheta1);
d2S_dtheta2_dtheta1 = 2*ddtheta2_dtheta2'*ddtheta2_ddtheta1;
d2S_dtheta2_2 = 2*(ddtheta2_dtheta2'*ddtheta2_dtheta2);

phi_idx = 1:nRoots;
dtheta1_idx = nRoots+1:nRoots+nDualEdges;
theta2_idx = nRoots+nDualEdges+1:nRoots+nDualEdges+nPixels;

H(phi_idx,phi_idx) =        wA*d2A_dphi2;
H(phi_idx,dtheta1_idx) =    wA*d2A_ddtheta1_dphi; 
H(dtheta1_idx, phi_idx) =   wA*d2A_ddtheta1_dphi';
H(phi_idx,theta2_idx) =     wA*d2A_dtheta2_dphi';
H(theta2_idx,phi_idx) =     wA*d2A_dtheta2_dphi; 
H(dtheta1_idx,dtheta1_idx) =wA*d2A_ddtheta1_2+wS*d2S_dtheta1_2;
H(dtheta1_idx,theta2_idx) = wA*d2A_ddtheta1_dtheta2' + wS*d2S_dtheta2_dtheta1'; 
H(theta2_idx,dtheta1_idx) = wA*d2A_ddtheta1_dtheta2 + wS*d2S_dtheta2_dtheta1; 
H(theta2_idx,theta2_idx) =  wA*d2A_dtheta2_2+wB*d2B_dtheta2_2+wS*d2S_dtheta2_2;

H = H + speye(size(H))*wDiff;
end

