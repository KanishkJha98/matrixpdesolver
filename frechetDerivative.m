function J = frechetDerivative(U,CH,CV,DH,DV,dt)
n = size(U,1);    
%Create the identity matrix of conforming dimensions
I = speye(n,n);
%Create a constant matrix that stores the value 1/dt
dT = (1/dt) * I;

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