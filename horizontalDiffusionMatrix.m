function [DH] = horizontalDiffusionMatrix(n,dy,Re)

% 5-point(3 order)central discretisation scheme for the discretisation
% of the diffusive term in the y-direction.The first and last rows are 0
% account for boundary nodes.
% function [DH] = verticalDiffusionMatrix_2Order(n,dy,Re)
% Input: n - The size of one dimension of the (square) matrix
%        dy - The size of the node spacing for y-axis, where dy = y(i+1) - y(i)
%        Re - The Reynold's number(reciprocal of the diffusion coefficient)
% Output: DH - A coefficient matrix that multiplies with U from the left to represent
%              the discretisation scheme for diffusive term in y-direction

i = [3:n-2 3:n-2 3:n-2 3:n-2 3:n-2];
j = [1:n-4 2:n-3 3:n-2 4:n-1 5:n];
v = [repelem((1/(12*dy^2))*(1/Re),n-4) repelem(-(16/(12*dy^2))*(1/Re),n-4) ...
    repelem((30/(12*dy^2))*(1/Re),n-4) repelem(-(16/(12*dy^2))*(1/Re),n-4)...
    repelem((1/(12*dy^2))*(1/Re),n-4)];
DH = sparse(i,j,v,n,n);

% 3-point(2 order) central discretisation scheme for the discretisation
% of nodes neighbouring to the boundary nodes
DH(2,[1,2,3]) = [-(1/dy^2)*(1/Re),(2/(dy^2))*(1/Re),-(1/dy^2)*(1/Re)];
DH(n-1,[n-2,n-1,n]) = [-(1/dy^2)*(1/Re),(2/(dy^2))*(1/Re),-(1/dy^2)*(1/Re)];

end