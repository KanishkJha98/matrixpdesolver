function [CH] = horizontalConvectiveMatrix(n,dy)

%5-point(4th order) central discretisation scheme for discretisation of
%convective term in y-direction.The top and bottom rows are 0
% to account for boundary nodes.
% function [CH] = horizontalConvectiveMatrix(n,dy)
% Input: n - The size of one dimension of the (square) matrix
%        dy - The size of the node spacing for y-axis, where dy = y(i+1) - y(i)
% Output: CH - A coefficient matrix that multiplies with U from the left to represent
%              the discretisation scheme for convective term in y-direction

i = [3:n-2 3:n-2 3:n-2 3:n-2];
j = [1:n-4 2:n-3 4:n-1 5:n];
v = [repelem((1/(12*dy)),n-4) repelem(-(8/(12*dy)),n-4) repelem((8/(12*dy)),n-4) ...
    repelem(-(1/(12*dy)),n-4)];
CH = sparse(i,j,v,n,n);
CH(2,[1,3]) = [-(1/(2*dy)),(1/(2*dy))];
CH(n-1,[n-2,n]) = [-(1/(2*dy)),(1/(2*dy))];

end