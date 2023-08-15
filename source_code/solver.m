clc

% Reynold's number, inverse of diffusion coefficient
R=100;
% Total time for which system is monitored
T = 0.25;
% The time step size
dt = 0.0005;

% The discretised mesh representing the region in which the system is
% studied
[X,Y] = meshgrid(0:1/49:1);
N = size(X,1);
h = 1/(N-1);

% The coefficient matrices on the nonlinear discretised system
CH = horizontalConvectiveMatrix_2Order(N,h);
CV = verticalConvectiveMatrix_2Order(N,h);
DH = horizontalDiffusionMatrix_2Order(N,h,R);
DV = verticalDiffusionMatrix_2Order(N,h,R);

%Assemble initial U as per initial value condition
Uold = 1./(1+exp((R*(X+Y)/2)));
Unew = Uold;
%Number of time steps
nt = T/dt;

%Current time
t = 0;

%tic and later toc are used to time the entire process, for gathering
%metrics
tic
% Compute eigendecomposition of DV to use for the Shifted Linear System
% solver for the Sylvester equation (COMMENT THE NEXT TWO CODE LINES FOR
%OTHER SOLVERS AS THIS WILL CONSUME TIME UNNECESSARILY)
[V,D] = eig(full(DV));
D = sparse(diag(D));

%Total number of linear iterations
lit = 0;
%Total number of Newton iterations
nit = 0;

% Loop through each time step and set the boundary values for the unknown U
% at the (k+1)th time step. Then compute the internal nodes for that time
% step from kth time step values using the Newton-Raphson method.
for k=1:nt
    t=t+dt;
    %Set boundary conditions at next time step for initial guess Unew
    Unew(1,1:N) = 1./(1+exp((R*(X(1,1:N)-t)/2)));
    Unew(2:N-1,1) = 1./(1+exp((R*(Y(2:N-1,1)-t)/2)));
    Unew(2:N-1,N) = 1./(1+exp((R*(Y(2:N-1,1)+1-t)/2)));
    Unew(N,1:N) = 1./(1+exp((R*(X(N,1:N)+1-t)/2)));
    fprintf("%f\n",t);
    [itlt,it,U,FU] = newtonSys( @(U)burgersDiscretisedForm_CN(U,Unew,Uold,CH,CV,DH,DV,dt), ...
         Unew, 1e-10, 30, ...
         @(U,FU)iterativeLinearSolveShiftedSylvester_CN(U,CH,CV,DH,DV,D,V,dt,FU,1e-10,100) );
    Unew = U;
    Uold = Unew;
    lit=lit+itlt;
    nit=nit+it;
end
totalTime = toc;

fprintf("Total number of linear iterations: %d\n",lit);
fprintf("Total number of Newton method iterations: %d\n",nit);
fprintf("Total time taken to run solver: %12.5f\n",totalTime);

% Uncomment to observe the flow described by the numerical solution to the
% equation
% figure(1)
% surf(X,Y,Unew)