function [J] = frechetDerivative_CN(U,CH,CV,DH,DV,dt)
% Standard approach to computing the Frechet derivative of the
% nonlinear matrix equation for the Crank-Nicolson scheme
%  function [J] = frechetDemoDerivative_CN( U,CH,CV,DH,DV,dt )
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
% Output: J - The Frechet Derivative equation, in vectorised form

n = size(U,1);    
%Create the identity matrix of conforming dimensions
I = speye(n,n);
%Create a constant matrix that stores the value 1/dt
dT = (2/dt) * I;

%Following operations can be done in parallel
M1 = kron(I,dT);
M2 = kron(I,DH);
M3 = kron(DV',I);
m4MatProd = sparse(CH*U);
M4 = diag(m4MatProd(:));
m5MatProd = sparse(U*CV);
M5 = diag(m5MatProd(:));

M6 = diag(sparse(U(:))) * kron(CV',I);
M7 = diag(sparse(U(:))) * kron(I,CH);
%Compute the Frechet derivative
J = M1 + M2 + M3 + M4 + M5 + M6 + M7;

end