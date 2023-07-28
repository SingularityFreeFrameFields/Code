function [b]=righthandside_b (numBoundaries, N, index_1, index_2)

%% Docstring
%
%    INPUT: 
%        numBoundaries: the total number of boundaries
%        N: the number of exterior boundaries
%        index_1/Exterior boundaries indices (EI): boundary tangent rotates 2pi, so index=(2pi-(angle sum))/2pi. It is a 1-by-(length(B)-N) matrix
%        index_2/Interior boundaries indices (II): boundary tangent rotates -2pi, so index=(-2pi-(angle sum))/2pi. It is a 1-by-(length(B)-N) matrix
%
%    OUTPUT:
%        b: a numBoundaries-by-1 vector containing the angle sum of all rotations the vector does per boundary      
%
%From the Poincare-Hopf theorem we have sum(all indices)=2-length(B). We should
%manually enter all indices in place of EIs and IIs in finalScript.m.

%Exterior boundaries 
p_1=[2*pi];
onevec_1= [ones(N,1)];
%index_1= [sym('II', [1 N])]'; this line is to show what the size of index_1 matrix is
b_1=(onevec_1-index_1)*p_1; 

%Interior boundaries
p_2=[2*pi];
onevec_2= [ones(numBoundaries-N, 1)];
%index_2= [sym('EI', [1 length(B)-N])]'; this line is to show what the size of index_2 matrix is
b_2=(-index_2-onevec_2)*p_2; 

b=[b_1;b_2];
