function J = frechetDemoDerivative(U,CH,CV,DH,DV,dt)
n = size(U,1);    
%Create the identity matrix of conforming dimensions
I = speye(n,n);
%Create a constant matrix that stores the value 1/dt
dT = (1/dt) * I;
%Following operations can be done in parallel
parPool = parpool("Threads",2);
f1 = parfeval(parPool,@kron,1,I,dT);
f2 = parfeval(parPool,@kron,1,I,DH);
f3 = parfeval(parPool,@kron,1,DV',I);
f4 = parfeval(parPool,@kron,1,CV',I);
f5 = parfeval(parPool,@kron,1,I,CH);

m4MatProd = sparse(CH*U);
M4 = gpuArray(diag(m4MatProd(:)));

m5MatProd = sparse(U*CV);
M5 = gpuArray(diag(m5MatProd(:)));

M1 = gpuArray(fetchOutputs(f1));
M2 = gpuArray(fetchOutputs(f2));
M3 = gpuArray(fetchOutputs(f3));
M6 = gpuArray(diag(sparse(U(:))) * fetchOutputs(f4));
M7 = gpuArray(diag(sparse(U(:))) * fetchOutputs(f5));
%Compute the Frechet derivative
J = M1 + M2 + M3 + M4 + M5 + M6 + M7;
delete(parPool);
end