function CCV = conservativeVerticalConvectiveMatrix(n)

% 3-point(2 order) central discretisation scheme for the discretisation
% of the conservative convective term in the x-direction. The first and
% last columns are 0 to account for boundary nodes.
% function [CCV] = conservativeHorizontalConvectiveMatrix(n,dy)
% Input: n - The size of one dimension of the (square) matrix
%        dx - The size of the node spacing for y-axis, where dy = y(i+1) - y(i)
% Output: CCV - A coefficient matrix that multiplies with U from the right to represent
%              the discretisation scheme for conservative convective term in x-direction

    i = [1:n-2 2:n-1];
    j = [2:n-1 2:n-1];
    v = [repelem((1/2),n-2) repelem((1/2),n-2)];
    CCV = sparse(i,j,v,n,n);
end