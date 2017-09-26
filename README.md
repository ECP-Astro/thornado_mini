# thornado_mini

This repository contains code to solve the equation of radiative transfer in the multi-group two-moment approximation.  
The Discontinuous Galekin (DG) method is used for spatial discretization, and an implicit-explicit (IMEX) method is used to integrate the moment equations in time.  The hyperbolic (streaming) part is treated explicitly, while the collision term is treated implicitly.

To compile and run the code you need a Fortran compiler, hdf5, mpi, and LAPACK.  

thornado_mini is divided into the main directories "Build", "Modules", "Workflow", and the application directory "DeleptonizationProblem1D".

The majority of the source code (i.e., solvers and associated utility functions) is under the directory "Modules". 

The directory "Workflow" contains the script "SetEnvironment.sh", which can be configured to set evnvironment variables, etc. for your machine, and Matlab scripts ("ReadFluidFields1D.m" and "ReadRadiationFields1D.m") that can be used to read simulation output for plotting.  
 
The application DeleptonizationProblem1D is described in Section 8.3 in ORNL/TM-2017/501 (Currently in review.)

"Build/Makefile_Compilers" contains machine specific compiler options, and is the only makefile that needs to be configured.  Example configurations are provided for compiling on a mac (using the gfortran compiler), and the titan machine at the Oak Ridge Leadership Computing Faccility (OLCF; using both the gnu and cray compilers).  Additional machines can be added by following the examples provided.  

Let's say you would like to compile DeleptonizationProblem1D on your machine using your favorite Fortran compiler.   We'll call this option mymachine:  
To do this, add lines in "Build/Makefile_Compilers" for FORTRAN_mymachine = , SUFFIX_f90_mymachine = , etc., and follow the instructions below.

To compile and run DeleptonizationProblem1D:

You are using the "mymachine" configuration discussed above, and you have placed the thornado_mini directory in "mydir"

1. Set the environment variable THORNADO_MACHINE:  

	export THORNADO_MACHINE=mymachine

2. Set the environment variable THORNADO_DIR:  

	export THORNADO_DIR=mydir/thornado_mini  

3. Go to mydir/thornado_mini/DeleptonizationProblem/Executables

	type "make".  This will generate the executable DeleptonizationProblem1D_mymachine

4. Copy or link the equation of state table (must be named named "EquationOfStateTable.h5") and the opacity table (must be named "OpacityTable.h5") into mydir/thornado_mini/DeleptonizationProblem/Executables.
   Equation of state tables and opacity tables needed for DeleptonizationProblem1D can be found at https://code.ornl.gov/m2o/thornado-mini-tables.  These are low-resolution tables: the equation of state table is about 53 MB, while the opacity table is about 260 MB.  

5. Run the executable in mydir/thornado_mini/DeleptonizationProblem/Executables

	./DeleptonizationProblem1D_mymachine

6. Output will be written to the directory mydir/thornado_mini/DeleptonizationProblem/Output

About the driver program Deleptonization1D:

The driver program is located in mydir/thornado_mini/DeleptonizationProblem/DeleptonizationProblem1D.f90.
It contains calls to InitializeProgram (to set up the mesh, allocate data structures, and initialize various parameters), InitializeDeleptonizationProblem1D (to set the initial condition for fluid and radiation fields), EvolveFields (to evolve the initial condition from t_begin to t_end, and write output after each time interval set by dt_write, using the Semi-Implicit Runge-Kutta (SIRK) method), and FinalizeProgram (to deallocate data structures).  

Notable parameters set in the call to InitializeProgram are:

nX_Option (INTEGER, DIMENSION(3)): Number of spatial elements in each spatial coordinate direction.  The code is currently implemented for one spatial dimension, so the second and third element should always be equal to one.  

xL_Option (REAL, DIMENSION(3)): Inner boundary coordinate in each spatial dimension.  

xR_Option (REAL, DIMENSION(3)): Outer boundary coordinate in each spatial dimension.  

zoomX_Option (REAL, DIMENSION(3)): Factor for geometrically increasing cell width, starting at the innermost element. A value of 1.0 results in an equidistant spatial grid in that dimension.  

nE_Option (INTEGER): Number of energy bins.  

eL_Option (REAL): Inner boundary of energy domain.  

eR_Option (REAL): Outer boundary of energy domain.  

zoomE_Option (REAL): Factor for geometrically increasing energy bin width, starting at the innermost bin. A value of 1.0 results in an equidistant energy grid.  

nNodes_Option (INTEGER): Number of nodes per active active dimension (i.e., with nX_Option > 1 or nE_Option > 1) in each element.  Polynomial degree = nNodes_Option - 1; i.e., nNodes_Option = 1 results in a spatially first order accurate method.  Possible values are nNodes_Option = 1, 2, 3, or 4.  

RadiationRiemannSolver_Option (CHARACTER): Riemann solver used to compute numberical fluxes across cell interfaces.  Available options are 'HLL' (the Harten-Lax-van Leer flux) and 'LLF' (the Local Lax-Friedrichs flux).  
