function [theta1_output] = fieldReconstruction(phi, d_theta1, theta2, narrowBand, edges_bfs, roots_bfs, pixelind_bfs, Tree, G, dualedgeInd)
[theta1, theta1_output] = integrate_theta(phi, d_theta1, edges_bfs,roots_bfs, pixelind_bfs, Tree, G, dualedgeInd, narrowBand);

hold on
axis equal

set(gca,'YDir','reverse');

q2 = quiver(Tree.Nodes.x(pixelind_bfs), Tree.Nodes.y(pixelind_bfs),cos(theta1/2),-sin(theta1/2),0.2,'g','linewidth',0.005);
q2.ShowArrowHead = 'off';
q3 = quiver(Tree.Nodes.x(pixelind_bfs), Tree.Nodes.y(pixelind_bfs),cos((theta2+theta1)/2),-sin((theta2+theta1)/2),0.2,'b','linewidth',0.005);    
q3.ShowArrowHead = 'off';
q4 = quiver(Tree.Nodes.x(pixelind_bfs), Tree.Nodes.y(pixelind_bfs),cos(theta1/2+pi),-sin(theta1/2+pi),0.2,'g','linewidth',0.005);    
q4.ShowArrowHead = 'off';
q5 = quiver(Tree.Nodes.x(pixelind_bfs), Tree.Nodes.y(pixelind_bfs),cos((theta2+theta1)/2+pi),-sin((theta2+theta1)/2+pi),0.2,'b','linewidth',0.005);
q5.ShowArrowHead = 'off';

end
