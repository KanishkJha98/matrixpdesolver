function [H] = iterativeLinearSolve( U,CH,CV,DH,DV,dt,FU,maxIt )

%Iterative GMRES solver for the Kronecker formulation approach of solving
%the linear matrix equation obtained from computing the Frechet derivative
% function [ H ] = iterativeLinearSolve( U,CH,CV,DH,DV,dt,FU,maxIt )
% Input: U - The matrix containing the unknown variables at every node
%            point
%        CH - The horizontal convective matrix corresponding to
%             discretisation of convective term in y-direction
%        CV - The vertical convective matrix corresponding to
%             discretisation of convective term in x-direction
%        DH - The horizontal diffusion matrix corresponding to
%             discretisation of the diffusive term in y-direction
%        DV - The vertical diffusion matrix corresponding to discretisation
%             of the diffusive term in x-direction
%        dt - The incremental time step value
%        FU - The R.H.S of the equation, which is just the nonlinear matrix
%        function that is to be finally computed
%        maxIt - maximum allowed number of iterations
% Output: H - The final solution of H or delta increment to U

n = size(U,1);
linTol = 1e-7;

FrD = frechetDerivative(U,CH,CV,DH,DV,dt);

% The chosen preconditioner is the solution to the nearest Kronecker
% product problem for the Frechet derivative, upto the 4th singular value
% iterate. Increasing the number of singular values used will lead to the
% preconditioner matrix being closer to the Frechet derivative, but also
% increase computational complexity drastically.
RA = reshape( permute( reshape( full(FrD), [ n, n, n, n ] ), [ 2 4 1 3 ] ), n^2, n^2 );
[URA,SRA,VRA] = svds(RA,4);
SRA = diag(SRA);
M = sparse(n^2,n^2);

parfor i=1:4
    M = M + sparse(SRA(i) .* kron(reshape(URA(:,i),n,n),reshape(VRA(:,i),n,n)));
end

if (isempty(M))
    [H,~,~,iter] = gmres(FrD,FU(:),[],linTol,maxIt);
else
    [H,~,~,iter] = gmres(FrD,FU(:),[],linTol,maxIt,M);
end

end