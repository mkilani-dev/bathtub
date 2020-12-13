# A code to solve for the optimum in the bathtub model

This project contains Fortran and AMPL code
to solve for the optimum in the Bathtub model.

Two types of utility functions are considered:

- alpha-beta-gamma tastes, and 
- smooth utility (logarithmic and exponential)

## Dependencies

The Fortran code calls some standard functions and libraries:

- DFZERO function
- lapack and Blas

and relies on 

- IpOpt (Fortran version)

for the optimization process.


