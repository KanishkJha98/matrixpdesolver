function [itlt,it, U,FU ] = newtonSys( fnon, Uinit, tol, maxIt, iterativeLinearSolve )

% Modified version of basic Newton-Raphson algorithm for solving a
% nonlinear matrix equation
%  function [ U,FU ] = newtonSys( fnon,Uinit, tol, maxit,iterativeLinearSolve )
% Input: fnon - function handle for nonlinear system
%        Uinit - Initial matrix of unknowns. Boundary nodes are set to
%        boundary values at the (k+1)th time step, i.e., the time step to
%        be computed.
%        tol - convergence tolerance
%        maxIt - maximum allowed number of iterations
%        iterativeLinearSolve - Function handle for an iterative linear
%        solver, to compute and solve the linear Frechet derivative equation obtained
% Output: U - Values at (k+1)th time step
%         FU - Function values of discretisation scheme at (k+1)th time
%         step. Boundary nodes are 0 as at boundary, FU = U(boundary)-U(boundary)
%         itlt - The total number of linear iterations conducted throughout
%                all the Newton iterations
%         it - The total number of Newton iterations

fprintf(' U      |f(U)|\n');

n = size(Uinit,1);
U = Uinit;
FU = feval(fnon,U);

it = 0;
itlt = 0;

fprintf(' %d %12.15f\n',it,norm(FU,"fro"));
%Compute Newton step till either norm of F(U) is less than tolerance, or we
%reach maximum number of allowed iterations
while (norm(FU,"fro")>tol && it < maxIt)

    %Solution of linear equation obtained from computing Frechet derivative
    [itl,H] = feval(iterativeLinearSolve,U,-FU);
    %Count number of linear iterations per Newton iteration and add them
    itlt = itlt + itl;
    %Update the values of H
    U = U + H;
    %Increment iteration counter by 1
    it = it + 1;
    
    FU = feval(fnon,U);
    fprintf(' %d %12.15f\n',it,norm(FU,"fro"));
end

if( it==maxIt)
    fprintf(' WARNING: Not converged\n')
else
    fprintf(' SUCCESS: Converged\n')
end
