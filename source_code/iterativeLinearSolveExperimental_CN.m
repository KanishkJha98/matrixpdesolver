function [it,H] = iterativeLinearSolveExperimental_CN( U,CH,CV,DH,DV,dt,FU,linTol,maxIt )

%An experimental solver algorithm which follows a splitting-based approach
% to solve the Frechet derivative equation,for the Crank-Nicolson scheme, 
% but does so by multiplying entire equation from the right by the kth 
% Identity vector.
%  function [ H ] = iterativeLinearSolveExperimental_CN( U,CH,CV,DH,DV,dt,FU,maxIt )
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
%        linTol - The allowed tolerance limit
%        maxIt - maximum allowed number of iterations
% Output: H - The final solution of H or delta increment to U 
%         it - The total number of linear iterations

n = size(U,1);

% Set initial guess to random values
H0 = rand(n,n);
H = zeros(n,n);

A1 = (2/dt) * speye(n,n);
A2 = CH * U;
A3 = U * CV;

it = 0;

RHSM = FU - H0 * DV - U .* (H0 * CV);
Resi = A1*H0 + H0 .* (CH * U) + U .* (CH * H0) + H0 .* (U * CV) + DH * H0 - RHSM;  

while (it < maxIt) && (norm(full(Resi))>linTol)
    parfor k=1:n
        LHS = A1 + sparse(diag(A2(:,k))) + sparse(diag(U(:,k))) * CH + sparse(diag(A3(:,k)))+ DH;
        H(:,k) = (LHS\RHSM(:,k));
    end
    it = it + 1;
    H0 = H;
    RHSM = FU - H0 * DV - U .* (H0 * CV);
    Resi = A1*H0 + H0 .* (CH * U) + U .* (CH * H0) + H0 .* (U * CV) + DH * H0 - RHSM;  
end
end