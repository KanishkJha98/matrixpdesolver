clc

% Reynold's number, inverse of diffusion coefficient
R=100;
% Total time for which system is monitored
T = 0.25;
% The time step size
dt = 0.0005;

% The discretised mesh representing the region in which the system is
% studied
[X,Y] = meshgrid(0:1/499:1);
N = size(X,1);
h = 1/(N-1);

% The coefficient matrices on the nonlinear discretised system
CH = horizontalConvectiveMatrix(N,h);
CV = verticalConvectiveMatrix(N,h);
DH = horizontalDiffusionMatrix(N,h,R);
DV = verticalDiffusionMatrix(N,h,R);
   

%Assemble initial U as per initial value condition
Uold = 1./(1+exp((R*(X+Y)/2)));

%Number of time steps
nt = T/dt;

%Current time
t = 0;

% Compute eigendecomposition of DV to use for the Shifted Linear System
% solver for the Sylvester equation
[V,D] = eig(full(DV));
D = sparse(diag(D));

% Loop through each time step and set the boundary values for the unknown U
% at the (k+1)th time step. Then compute the internal nodes for that time
% step from kth time step values using the Newton-Raphson method.
for k=1:nt
    t=t+dt;
    %Set boundary conditions at every time step
    Uold(1,1:N) = 1./(1+exp((R*(X(1,1:N)-t)/2)));
    Uold(2:N-1,1:N-1:N) = 1./(1+exp((R*(X(2:N-1,1:N-1:N)+Y(2:N-1,1:N-1:N)-t)/2)));
    Uold(N,1:N) = 1./(1+exp((R*(X(N,1:N)+1-t)/2)));
    fprintf("%f\n",t);
    tic
    [U,H,FU] = newtonSys( @(U)burgersDiscretisedForm(U,Uold,CH,CV,DH,DV,dt), ...
         Uold, 1e-12, 20, ...
         @(U,FU)iterativeLinearSolveShiftedSylvester(U,CH,CV,DH,D,V,dt,FU,100) );
    toc
    Uold = U;

end


figure(1)
surf(X,Y,Uold)

temp = 1./(1+exp((R*(X+Y-t)/2)));

figure(4)
surf(X,Y,temp)