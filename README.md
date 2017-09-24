# thornado_mini

This repository contains code to solve the equation of radiative transfer in the multi-group two-moment approximation.  
The Discontinuous Galekin (DG) method is used for spatial discretization, and an implicit-explicit (IMEX) method is used to integrate the moment equations in time.  The hyperbolic (streaming) part is treated explicitly, while the collision term is treated implicitly.

To compile and run the code you need a Fortran compiler, hdf5, mpi, and LAPACK.  

thornado_mini is divided into the main directories "Build", "Modules", "Workflow", and the application directory "DeleptonizationProblem1D".

"Build/Makefile_Compilers" contains machine specific compiler options, and is the only makefile that needs to be configured.  Currently, configurations are provided as examples for compiling on a mac (using the gfortran compiler), and the titan machine at the Oak Ridge Leadership Computing Faccility (OLCF; using both the gnu and cray compilers).  Additional machines can be added by following the templates provided.

The majority of the source code (i.e., solvers and associated utility functions) is under the directory "Modules".

"Workflow" contains the script "SetEnvironment.sh" and Matlab scripts ("ReadFluidFields1D.m" and "ReadRadiationFields1D.m") that can be used to read simulation output for plotting.  

How to compile and run DeleptonizationProblem1D:
