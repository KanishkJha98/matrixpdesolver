function [H] = iterativeLinearSolveShiftedSylvester( U,CH,CV,DH,DV,D,V,dt,FU,linTol,maxIt )

%Iterative split-method relaxation-like solver for the Frechet derivative
%using Shifted linear system solve approach for resulting Sylvester
%equation
%  function [ H ] = iterativeLinearSolveShiftedSylvester( U,CH,CV,DH,DV,D,V,dt,FU,linTol,maxIt )
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
% Set initial guess to random values
H0 = rand(n,n);
H = zeros(n,n);
A1 = (1/dt) * speye(n,n);
A2 = sparse(CH * U);
A3 = sparse(U * CV);
it=0;
RHS = FU - H0 .* A2 - U .* (CH * H0) - H0 * A3 - U .* (H0*CV);
Resi = (A1+DH)*H + H*DV - RHS;
while (it < maxIt && norm(Resi,"fro") > linTol)
    it=it+1;
    H = sylvesterSolver(A1+DH,D,V,RHS);
    H0 = H;
    RHS = FU - H0 .* A2 - U .* (CH * H0) - H0 * A3 - U .* (H0*CV);
    Resi = (A1+DH)*H + H*DV - RHS;
end
fprintf("%d\n",it);
end