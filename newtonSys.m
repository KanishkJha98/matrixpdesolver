function [ U,FU ] = newtonSys( fnon, Uold, tol, maxIt, iterativeLinearSolve )

% Modified version of basic Newton-Raphson algorithm for solving a
% nonlinear matrix equation
%  function [ U,FU ] = newtonSys( fnon,Uold, tol, maxit,iterativeLinearSolve )
% Input: fnon - function handle for nonlinear system
%        Uold - Initial matrix of unknowns. Boundary nodes are set to
%        boundary values at the (k+1)th time step, i.e., the time step to
%        be computed.
%        tol - convergence tolerance
%        maxIt - maximum allowed number of iterations
%        iterativeLinearSolve - Function handle for an iterative linear
%        solver, to compute and solve the linear Frechet derivative equation obtained
% Output: U - Values at (k+1)th time step
%         FU - Function values of discretisation scheme at (k+1)th time
%         step. Boundary nodes are 0 as at boundary, FU = U(boundary)-U(boundary)

fprintf(' U      |f(U)|\n');

n = size(Uold,1);
U = Uold;            
FU = feval(fnon,U);
it = 0;

fprintf(' %d %12.15f\n',it,norm(FU,"fro"));
%Compute Newton step till either norm of F(U) is less than tolerance, or we
%reach maximum number of allowed iterations
while (norm(FU,"fro")>tol && it < maxIt)
    %Solution of linear equation obtained from computing Frechet derivative
    [H] = feval(iterativeLinearSolve,U,-FU);
    %Set boundary values to 0, as we do not wish to change the U values at
    %the boundary nodes.
    H(1,1:n) = 0;
    H(2:n-1,1:n-1:n) = 0;
    H(n,1:n) = 0;
    %Update the values of H
    U = U + H;
    FU = feval(fnon,U);

    it = it + 1;
    fprintf(' %d %12.15f\n',it,norm(FU,"fro"));
end

if( it==maxIt)
    fprintf(' WARNING: Not converged\n')
else
    fprintf(' SUCCESS: Converged\n')
end
