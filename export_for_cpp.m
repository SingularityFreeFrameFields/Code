function []= export_for_cpp(filename,phi,d_theta1,theta2,wA, wB, wS, antialign, ddtheta2_dtheta2, dtheta1_ddtheta1, ddtheta2_dtheta1, edges_bfs, roots_bfs, pixelind_bfs, ind_of_zero_gradients_BFS, Tree, dualEdgeInd, all_edges_bfs, myBFSorder, G, PixelIndex, g_theta_BFS, finalA, finalb)
f = fopen(filename,'w');
name("phi");
saveVector(phi);
name("d_theta1");
saveVector(d_theta1);
name("theta2");
saveVector(theta2);
name("wA");
fprintf(f,"%.9f\n",wA);
name("wB");
fprintf(f,"%.9f\n",wB);
name("wS");
fprintf(f,"%.9f\n",wS);
name("antialign");
fprintf(f,"%.9f\n",antialign);

name("ddtheta2_dtheta2");
saveSparseMatrix(ddtheta2_dtheta2);
name("dtheta1_ddtheta1");
saveSparseMatrix(dtheta1_ddtheta1);
name("ddtheta2_dtheta1");
saveSparseMatrix(ddtheta2_dtheta1);

name("edges_bfs");
savePairs(edges_bfs);
name("roots_bfs");
saveVector(int32(roots_bfs));
name("pixelind_bfs");
saveVector(int32(pixelind_bfs));
name("ind_of_zero_gradients_BFS");
saveVector(int32(ind_of_zero_gradients_BFS));

name("Tree");
savePairs(Tree.Edges.EndNodes);
name("dualEdgeInd");
savePairs(dualEdgeInd);
name("alledges_bfs");
savePairs(all_edges_bfs);

name("myBFSorder");
saveVector(int32(myBFSorder));
name("G");
savePairs(G.Edges.EndNodes);
name("PixelIndex");
saveVector(int32(PixelIndex));
name("g_theta_BFS");
saveVector(g_theta_BFS);

name("finalA");
saveSparseMatrix(finalA);
name("finalb");
saveVector(finalb);

fclose(f);

    function [] = saveVector(v)
        fprintf(f, "%d\n", numel(v));
        for jj=1:numel(v)
            if isa(v(jj),'double')
                fprintf(f,"%.9f\n",v(jj));
            else
                fprintf(f,"%d\n",v(jj)-1); %for indices to start from 0
            end
        end
    end

    function [] = saveSparseMatrix(M)
     fprintf(f,"%d %d %d\n", nnz(M), size(M,1), size(M,2));
     [ii,jj,val] = find(M);
     for k = 1:numel(ii)
         fprintf(f,"%d %d %.9f\n",ii(k)-1,jj(k)-1,val(k));
     end
    end

    function [] = savePairs(p)
        fprintf(f, "%d\n", size(p,1));
        for jj=1:size(p,1)
            fprintf(f,"%d %d\n",p(jj,1)-1,p(jj,2)-1);
        end
    end

    function [] = name(str)
      fprintf(f,"%s\n",str);
    end
end