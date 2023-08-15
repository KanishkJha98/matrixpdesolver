function [DV] = verticalDiffusionMatrix_2Order(n,dx,Re)

% 3-point(or 2nd Order) central discretisation scheme for discretisation
% of the diffusive term in the x-direction. The first and last columns are 0
% account for boundary nodes. 
% function [DV] = verticalDiffusionMatrix_2Order(n,dx,Re)
% Input: n - The size of one dimension of the (square) matrix
%        dx - The size of the node spacing for x-axis, where dx = x(i+1) - x(i)
%        Re - The Reynold's number(reciprocal of the diffusion coefficient)
% Output: DV - A coefficient matrix that multiplies with U from the right to represent
%              the discretisation scheme for diffusive term in x-direction

i = [1:n-2 2:n-1 3:n];
j = [2:n-1 2:n-1 2:n-1];
v = [repelem(-(1/(dx^2))*(1/Re),n-2) repelem((1/(dx^2))*(2/Re),n-2) repelem(-(1/(dx^2))*(1/Re),n-2)];
DV = sparse(i,j,v,n,n);
end