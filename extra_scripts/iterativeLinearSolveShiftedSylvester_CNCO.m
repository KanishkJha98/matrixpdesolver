function [H] = iterativeLinearSolveShiftedSylvester_CNCO( U,CH,CCH,CV,CCV,DH,DV,D,V,dt,FU,linTol,maxIt )

%Iterative split-method relaxation-like solver for the Frechet derivative
%using Shifted linear system solve approach for resulting Sylvester
%equation, for the Crank-Nicolson scheme
%  function [ H ] = iterativeLinearSolveShiftedSylvester_CN( U,CH,CV,DH,DV,D,V,dt,FU,linTol,maxIt )
% Input: U - The matrix containing the unknown variables at every node
%            point
%        CH - The horizontal convective matrix corresponding to
%             discretisation of convective term in y-direction
%        CCH - A coefficient matrix that multiplies with U from the left to represent
%              the discretisation scheme for conservative convective term in y-direction
%        CV - The vertical convective matrix corresponding to
%             discretisation of convective term in x-direction
%        CCV - A coefficient matrix that multiplies with U from the right to represent
%              the discretisation scheme for conservative convective term in x-direction
%        DH - The horizontal diffusion matrix corresponding to
%             discretisation of the diffusive term in y-direction
%        DV - The vertical diffusion matrix corresponding to discretisation
%             of the diffusive term in x-direction
%        D - The diagonal vector containing eigenvalues of the vertical
%            diffusion matrix corresponding to discretisation of the
%            diffusive term in x-direction
%        V - The eigenvector matrix corresponding to eigendecomposition of
%            the vertical diffusion matrix
%        dt - The incremental time step value
%        FU - The R.H.S of the equation, which is just the nonlinear matrix
%        function that is to be finally computed
%        linTol - The accepted tolerance limit
%        maxIt - maximum allowed number of iterations
% Output: H - The final solution of H or delta increment to U

n = size(U,1);
%Ensure H and H0 are not initially 0, as this will lead to the convergence
%check being successful instantly.
H = zeros(n,n);
H0 = zeros(n,n);
A1 = (2/dt) * speye(n,n);
A2 = sparse(CH * U);
A3 = sparse(U * CV);
it=0;
RHS = FU - (CCH*H0) .* A2 - (CCH*U) .* (CH * H0) - (H0*CCV) .* A3 - (U*CCV) .* (H0*CV);
Resi = (A1+DH)*H + H*DV - RHS;
while (it < maxIt && norm(Resi,"fro") > linTol)
    it=it+1;
    H = sylvesterSolver(A1+DH,D,V,RHS);
    H0 = H;
    RHS = FU - (CCH*H0) .* A2 - (CCH*U) .* (CH * H0) - (H0*CCV) .* A3 - (U*CCV) .* (H0*CV);
    Resi = (A1+DH)*H + H*DV - RHS;
end
fprintf("%d\n",it);
end