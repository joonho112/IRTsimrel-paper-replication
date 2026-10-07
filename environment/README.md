<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# Software environment

This release was exercised with R 4.6.0 on macOS. `packages.csv` and `renv.lock` record the dependency versions available during the release checks. These files describe the checked environment; they do not imply that every historical job used exactly those versions.

The retained fixed-form and node-study computations used IRTsimrel 0.3.0. The worked example and item-population computations used 0.3.1. `code/install_dependencies.R` installs both source archives into `.library/0.3.0/` and `.library/0.3.1/`. Every study verifies the requested namespace and refuses to use a different loaded version. Run different versions in separate R processes.

The convenience installer installs missing CRAN dependencies. It does not downgrade existing packages. To restore the recorded CRAN dependency versions into a separate project library, install `renv` and run:

```r
renv::restore(lockfile = "renv.lock", library = ".library/cran", prompt = FALSE)
```

Then use `R_LIBS_USER=.library/cran` when launching R, and install the two frozen archives with `Rscript code/install_dependencies.R --packages-only`. A fresh cross-platform restore was not part of the release checks. Compiled dependencies such as `lme4` and `TAM` may require R development tools when binaries are unavailable.

We pin `Mersenne-Twister`, `Inversion`, and `Rejection` for the RNG. Numerical workers request one BLAS/OpenMP thread each. Some BLAS implementations read thread settings only at process startup; the shell launch environment can also set `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` before running R.

Exact reproducibility is strongest for the calibration kernels. Mixed-model fits can vary near boundary optima across platforms and library versions. Inspect statuses and discrepancies; do not remove them to force agreement.
