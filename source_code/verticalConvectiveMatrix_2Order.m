function [CV] = verticalConvectiveMatrix_2Order(n,dx)

% 3-point(or 2nd order) accurate central discretisation scheme for
% discretisation of convective term in x-direction. The first and last
% columns are 0 to account for boundary nodes.
% function [CV] = verticalConvectiveMatrix_2Order(n,dx)
% Input: n - The size of one dimension of the (square) matrix
%        dx - The size of the node spacing for x-axis, where dx = x(i+1) - x(i)
% Output: CV - A coefficient matrix that multiplies with U from the right to represent
%              the discretisation scheme for convective term in x-direction

j = [2:n-1 2:n-1];
i = [1:n-2 3:n];
v = [repelem(-(1/(2*dx)),n-2) repelem((1/(2*dx)),n-2)];
CV = sparse(i,j,v,n,n);

end