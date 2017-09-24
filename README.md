# thornado_mini

This repository contains code to solve the equation of radiative transfer in the multi-group two-moment approximation.  
The Discontinuous Galekin (DG) method is used for spatial discretization, and an implicit-explicit (IMEX) method is used to integrate the moment equations in time.  The hyperbolic (streaming) part is treated explicitly, while the collision term is treated implicitly.

To compile and run the code you need a Fortran compiler, hdf5, mpi, and LAPACK.  

thornado_mini is divided into the main directories "Build", "Modules", "Workflow", and the application directory "DeleptonizationProblem1D".

DeleptonizationProblem1D is described in Section 8.3 in ORNL/TM-2017/501 (Currently in review.)

"Build/Makefile_Compilers" contains machine specific compiler options, and is the only makefile that needs to be configured.  Example configurations are provided for compiling on a mac (using the gfortran compiler), and the titan machine at the Oak Ridge Leadership Computing Faccility (OLCF; using both the gnu and cray compilers).  Additional machines can be added by following the examples provided.  

Let's say you would like to compile using your favorite compiler.   We'll call this option mymachine:  
Then add lines in "Build/Makefile_Compilers" for FORTRAN_mymachine = , SUFFIX_f90_mymachine = , etc.

The majority of the source code (i.e., solvers and associated utility functions) is under the directory "Modules".

The directory "Workflow" contains the script "SetEnvironment.sh", which can be configured to set evnvironment variables, etc., and Matlab scripts ("ReadFluidFields1D.m" and "ReadRadiationFields1D.m") that can be used to read simulation output for plotting.  

To compile and run DeleptonizationProblem1D:

You are using the "mymachine" configuration discussed above, and you have placed the thornado_mini directory in "mydir"

1. Set the environment variable THORNADO_MACHINE:  

	export THORNADO_MACHINE=mymachine

2. Set the environment variable THORNADO_DIR:  

	export THORNADO_DIR=mydir/thornado_mini  

3. Go to mydir/thornado_mini/DeleptonizationProblem/Executables

	type "make".  This will generate the executable DeleptonizationProblem1D_mymachine

4. Copy or link the equation of state table (must be named named "EquationOfStateTable.h5") and the opacity table (must be named "OpacityTable.h5") into mydir/thornado_mini/DeleptonizationProblem/Executables

5. Run the executable in mydir/thornado_mini/DeleptonizationProblem/Executables

	./DeleptonizationProblem1D_mymachine

6. Output will be written to the directory mydir/thornado_mini/DeleptonizationProblem/Output
