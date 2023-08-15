function [DH] = horizontalDiffusionMatrix_2Order(n,dy,Re)

%3-point(2 order) central discretisation scheme for discretisation of the
%diffusion term in the y-direction.The first and last rows are 0
% account for boundary nodes. 
% function [DH] = verticalDiffusionMatrix_2Order(n,dy,Re)
% Input: n - The size of one dimension of the (square) matrix
%        dy - The size of the node spacing for y-axis, where dy = y(i+1) - y(i)
%        Re - The Reynold's number(reciprocal of the diffusion coefficient)
% Output: DH - A coefficient matrix that multiplies with U from the left to represent
%              the discretisation scheme for diffusive term in y-direction

i = [2:n-1 2:n-1 2:n-1];
j = [1:n-2 2:n-1 3:n];
v = [repelem(-(1/(dy^2))*(1/Re),n-2) repelem((1/(dy^2))*(2/Re),n-2) repelem(-(1/(dy^2))*(1/Re),n-2)];
DH = sparse(i,j,v,n,n);

end