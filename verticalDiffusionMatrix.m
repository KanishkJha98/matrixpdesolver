function [DV] = verticalDiffusionMatrix(n,dx,Re)

% 5-point(or 4th Order) central discretisation scheme for discretisation
% of the diffusive term in the x-direction. The first and last columns are 0
% to account for boundary nodes.
% function [DV] = verticalDiffusionMatrix(n,dx,Re)
% Input: n - The size of one dimension of the (square) matrix
%        dx - The size of the node spacing for x-axis, where dx = x(i+1) - x(i)
%        Re - The Reynold's number(reciprocal of the diffusion coefficient)
% Output: DV - A coefficient matrix that multiplies with U from the right to represent
%              the discretisation scheme for diffusive term in x-direction

j = [3:n-2 3:n-2 3:n-2 3:n-2 3:n-2];
i = [1:n-4 2:n-3 3:n-2 4:n-1 5:n];
v = [repelem((1/(12*dx^2))*(1/Re),n-4) repelem(-(16/(12*dx^2))*(1/Re),n-4) ...
     repelem((30/(12*dx^2))*(1/Re),n-4) repelem(-(16/(12*dx^2))*(1/Re),n-4)...
     repelem((1/(12*dx^2))*(1/Re),n-4)];
DV = sparse(i,j,v,n,n);
    
% 3-point(or 2nd Order) central discretisation scheme for discretisation
% of nodes neighbouring to the boundary nodes.
DV([1,2,3],2) = [-(1/dx^2)*(1/Re),(2/(dx^2))*(1/Re),-(1/dx^2)*(1/Re)];
DV([n-2,n-1,n],n-1) = [-(1/dx^2)*(1/Re),(2/(dx^2))*(1/Re),-(1/dx^2)*(1/Re)];

end