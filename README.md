# thornado_mini

This repository contains code to solve the equation of radiative transfer in the two-moment approximation.  
The Discontinuous Galekin (DG) method is used for spatial discretization, and an implicit-explicit (IMEX) method is used to integrate the moment equations in time.  The hyperbolic (streaming) part is treated explicitly, while the collision term is treated implicitly.  