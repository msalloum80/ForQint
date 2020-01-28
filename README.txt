This folder contains the fortran quadrature library and a test case.
====================================================================

- Build the library and the executable by typing "make" in the source directory

- run "run.x" to see example output

What this library contains
=============================

lib - directory which stores the fortran routine to generate quadrature points.

quad.F90 - contains routines to compute integrals via quadrature (this is the main ForQint routines).

fnfm.F90 - function implementations of updraft velocity, i.e., function implementation to calculate fn, fm, fluxfn and fluxfm (for test case).

main.F90 - test code to integrate functions defined in fnfm. Note that some parameters (e.g. alpha, beta...) were computed beforehand.