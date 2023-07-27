function [CV] = verticalConvectiveMatrix(n,dx)

% 5-point(or 3rd order accurate) central discretisation scheme for the
% discretisation of the convective term in the x-direction. The first and last
% columns are 0 to account for boundary nodes.
% function [CV] = verticalConvectiveMatrix(n,dx)
% Input: n - The size of one dimension of the (square) matrix
%        dx - The size of the node spacing for x-axis, where dx = x(i+1) - x(i)
% Output: CV - A coefficient matrix that multiplies with U from the right to represent
%              the discretisation scheme for convective term in x-direction

i = [1:n-4 2:n-3 4:n-1 5:n];
j = [3:n-2 3:n-2 3:n-2 3:n-2];
v = [repelem((1/(12*dx)),n-4) repelem(-(8/(12*dx)),n-4) repelem((8/(12*dx)),n-4) ...
     repelem(-(1/(12*dx)),n-4)];
CV = sparse(i,j,v,n,n);
    
% For nodes neighbouring the boundary nodes, apply 3-point or 2nd order
% accurate central finite difference scheme
CV([1,3],2) = [-(1/(2*dx)),(1/(2*dx))];
CV([n-2,n],n-1) = [-(1/(2*dx)),(1/(2*dx))];

end