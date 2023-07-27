function [FU] = burgersDiscretisedForm(U,Uold,CH,CV,DH,DV,dt)

% Store the discretised form of the Burgers' equation as a nonlinear matrix
% equations. The boundary values are set to 0 after we use them in the equation,
% as we do not necessarily need the data at the boundary nodes for any form
% of computation. They exist merely to ensure that the dimensions of the
% coefficient matrices conform.
% function [FU] = burgersDiscretisedForm(U,Uold,CH,CV,DH,DV,dt)
% Input: U - The matrix containing the unknown variables at every node
%            point
%        Uold - The matrix containing the (known) values of U at current
%        time step
%        CH - The horizontal convective matrix corresponding to
%             discretisation of convective term in y-direction
%        CV - The vertical convective matrix corresponding to
%             discretisation of convective term in x-direction
%        DH - The horizontal diffusion matrix corresponding to
%             discretisation of the diffusive term in y-direction
%        DV - The vertical diffusion matrix corresponding to discretisation
%             of the diffusive term in x-direction
%        dt - The incremental time step value
% Output: FU - The equation matrix that represents the discretised form of
%              the Burgers' equation for each internal node, and has 0 at
%              each boundary node.

n = size(U,1);
FU = (1/dt) * U - (1/dt) * Uold + U .* (U * CV) + U .* (CH * U) + U * DV + DH * U;

%Set boundary values of FU equal to 0(can
%parallelise)
FU(1,1:n) = U(1,1:n)-Uold(1,1:n);
FU(2:n-1,1:n-1:n) = U(2:n-1,1:n-1:n)-Uold(2:n-1,1:n-1:n);
FU(n,1:n) = U(n,1:n)-Uold(n,1:n);

end