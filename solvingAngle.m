function [solangles,finalA,finalb]=SolvingAngle(H, d0sub, b, singularInit)

finalA=[[d0sub]';H] ;

b_1=zeros(size(d0sub,2),1);

finalb=[b_1;b];
min_norm_from_singular = lsqminnorm(finalA, finalb - finalA*singularInit);

solangles = singularInit + min_norm_from_singular;
