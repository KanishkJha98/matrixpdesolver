function [H] = iterativeLinearSolveMatlabSylvester_CNCO( U,CH,CCH,CV,CCV,DH,DV,dt, FU,linTol,maxIt )

%Iterative split-method relaxation-like solver for the Frechet derivative
%using Matlab's sylvester function to solve resulting Sylvester
%equation, for the Crank-Nicolson scheme.
% function [ H ] = iterativeLinearSolveShiftedSylvester_CN( U,CH,CV,DH,DV,dt,FU,maxIt )
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
%        dt - The incremental time step value
%        FU - The R.H.S of the equation, which is just the nonlinear matrix
%        function that is to be finally computed
%        linTol - The accepted tolerance limit
%        maxIt - maximum allowed number of iterations
% Output: H - The final solution of H or delta increment to U

n = size(U,1);
%Initialise the initial guess to a random set of variables. 
H0 = rand(n,n);
H = zeros(n,n);
A1 = (2/dt) * speye(n,n);
A2 = sparse(CH * U);
A3 = sparse(U * CV);
RHS = FU - (CCH*H0) .* A2 - (CCH*U) .* (CH * H0) - (H0*CCV) .* A3 - (U*CCV) .* (H0*CV);
Resi = (A1+DH)*H + H*DV - RHS;
it=0;
while (it < maxIt && norm(Resi,"fro") > linTol)
    it=it+1;
    H = sylvester(full(A1+DH),full(DV),RHS);
    H0 = H;
    RHS = FU - (CCH*H0) .* A2 - (CCH*U) .* (CH * H0) - (H0*CCV) .* A3 - (U*CCV) .* (H0*CV);
    Resi = (A1+DH)*H + H*DV - RHS;
end

end