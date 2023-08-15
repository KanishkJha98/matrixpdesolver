The code files in this repository are designed to produce a blueprint and a starting point towards solving complex 2-dimensional nonlinear Partial Differential Equations via simple Finite Difference Schemes, but without vectorising the matrix of unknowns. Instead, through careful manipulation and formulation of equivalent equations, the codebase is geared towards utilising the matrix structure of the problem in solving the equations.<br> Similar approaches can be found with regards to linear elliptic PDEs, advection-diffusion equations and parabolic PDEs. However, as far as my knowledge and research goes, approaches geared towards solving nonlinear PDEs are nearly nonexistent. Moreover, the few research articles that do exist result in the formulation of nonlinear matrix equations that only exist in two well-suited operators: the standard matrix product and the standard matrix addition. This particular problem, that stems from a variant of the 2-dimensional viscous Burgers' equation, results in a matrix equation that, on top of being nonlinear, contains both Hadamard and matrix products, alongside matrix addition. This codebase is my attempt at providing certain solvers to still utilise matrix structure in these scenarios to solve the problem.<br> Note that much improvement may and can be done to the code. This is merely a first attempt at understanding this problem, and ensuring that this particular problem is brought forward in a cohesive manner. The equation being solved is $$\frac {\partial u}{\partial t} + u \frac {\partial u}{\partial x} + u \frac {\partial u}{\partial y} = \nu (\frac {\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2})$$ The corresponding nonlinear matrix equation is $$\frac{2}{\Delta t} . U - \frac{2}{\Delta t} . Uold + DH . U + (CH.U) \odot U + U . DV + (U . CV) \odot V + DH . Uold + (CH.Uold) \odot Uold + Uold . DV + (Uold . CV) \odot V = 0$$ where the space discretisation scheme used is a 3-point central finite difference scheme, and the discretisation scheme in time is the Crank-Nicolson method. U is the matrix of unknown variables over the entire grid being studied, and all other matrices represent the coefficients of discretisation. Here $\odot$ is used to represent the Hadamard product, and $.$ to represent the normal matrix product.<br>
There are two separate folders containing script files. The first folder titled "source_code" contains all the scripts used for the aforementioned schemes and general solvers. This folder represents the code used in Chapters 4 and 5 of the corresponding report. The solver.m script file is used to run all the other script files, and contains the initial conditions, different parameter values like Reynold's number, time-step size and grid size. It also performs runtime measurement. It is the starting point of the entire solver, as it contains the loop for time-stepping, within which boundary conditions are set at each time step followed by application of the Newton method. The second folder title "extra_scripts" contain code mentioned in the Appendix of the report. This includes alternative 4th order Central Finite Difference Scheme, Implicit Euler scheme, and the necessary matrices for solving the Conservative form of the equation.
