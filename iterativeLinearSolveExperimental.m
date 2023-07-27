function [H] = iterativeLinearSolveExperimental( U,CH,CV,DH,DV,dt,FU,maxIt )

%An experimental solver algorithm which follows a splitting-based approach
% to solve the Frechet derivative equation but does so by multiplying 
% entire equation from the right by the kth Identity vector. Convergence 
% for the method depends on utilising some learning parameter omega.
%  function [ H ] = iterativeLinearSolveExperimental( U,CH,CV,DH,DV,dt,FU,maxIt )
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
omega = 0.5;

H = zeros(n,n);
H0 = zeros(n,n);

A1 = (1/dt) * speye(n,n);
A2 = sparse(CH * U);
A3 = sparse(U * CV);

RHSM = zeros(n,n);
it = 0;
parfor k=1:n
    RHSM(:,k) = FU(:,k) - H0 * DV(:,k) - sparse(diag(U(:,k))) * H0 * CV(:,k);
end
while (it < maxIt) && (norm(full(RHSM))>linTol)
    parfor k=1:n
        LHS = A1 + sparse(diag(A2(:,k))) + sparse(diag(U(:,k))) * CH + sparse(diag(A3(:,k)))+ DH;
        
        H(:,k) = omega*(LHS\RHSM(:,k));
        
        it = it + 1;
    end
    H0 = H;
    RHSM = FU - H0 * DV - U .* (H0 * CV);
end

end