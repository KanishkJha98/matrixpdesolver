function [X] = sylvesterSolver(A,D,V,C)

% Solves the Sylvester equation via solving a system of shifted linear
% systems for every column of a permuted solution matrix X*, and then
% computes the actual solution matrix by multiplying it with the inverse of
% the permutation matrix or eigenvector matrix V.
% function [X] = sylvesterSolver(A,D,V,C)
% Input: A - The first coefficient matrix, that multiplies X from the left
%        D - A column matrix containing the eigenvalues of the second
%        coefficient matrix B that multiplies X from the right
%        V - The eigenvector matrix of B
%        C - The matrix of constant values, that exists as the R.H.S of the
%        Sylvester equation
% Output: X - Solution matrix of the Sylvester equation

n = size(A,1);
X = zeros(n,n);
I = speye(n,n);
C = C*V;
parfor k=1:n
    X(:,k) = (A+D(k)*I)\C(:,k);
end
X = X/V;
end