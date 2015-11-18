Approximations For Solving L0-Norm Minimization Problems In Compressed Sensing
==============================================================================

Compressed Sensing is a signal processing technique for under-sampling and optimal recovery of a signal, in particular images. It exploits the fact that natural images are sparse in some domain like the Discrete Fourier Transform or Wavelets. The sparsity of the signal can be exploited via optimization to recover it from far fewer samples than required by the Nyquist-Shannon sampling theorem. The sparsity constraint, however, is of a combinatorial nature and is not convex. We explore various methods for solving this problem via convex optimization.

Running the code
================

-   OMP: Go to the `omp/` folder and set the `m` parameter in `omp.m` to specify fraction of available measurements. Then run `omp.m`.
-   SCA: Go to the `sca/` folder and set the `m` parameter in `sca.m` or `scap.m` to specify fraction of available measurements, and the `k`, `t`, `eps`, and `maxIter` parameters for controlling sparsity of the output, coarseness of the L0 norm approximation, stopping error and maximum iterations per frame respectively. Then run `sca.m` or `scap.m`. `scap.m` is the parallel version of `sca.m`.
-   PD: Go to the `pd/` folder and set the `m` parameter in `pd_parallel.m` to specify fraction of available measurements, and other parameters to run the PD reconstruction in the solvePD.m code itself. 

All code has been written and tested on Matlab R2014b.

Authors
=======
Alankar Kotwal, <alankarkotwal13@gmail.com>  
Anand Kalvit, <anandiitb12@gmail.com>
