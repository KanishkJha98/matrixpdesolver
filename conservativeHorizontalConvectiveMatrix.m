function [CCH] = conservativeHorizontalConvectiveMatrix(n)

% 3-point(2 order) central discretisation scheme for the discretisation
% of the conservative convective term in the y-direction. The top and 
% bottom rows are 0 to account for boundary nodes.
% function [CCH] = conservativeHorizontalConvectiveMatrix(n,dy)
% Input: n - The size of one dimension of the (square) matrix
%        dy - The size of the node spacing for y-axis, where dy = y(i+1) - y(i)
% Output: CCH - A coefficient matrix that multiplies with U from the left to represent
%              the discretisation scheme for conservative convective term in y-direction

    i = [2:n-1 2:n-1];
    j = [1:n-2 2:n-1];
    v = [repelem((1/2),n-2) repelem((1/2),n-2)];
    CCH = sparse(i,j,v,n,n);
end