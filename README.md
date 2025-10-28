# ContinuumSystem

Code used to compute the mean-field data and collective excitation spectra presented in 

Collective modes in superconductors including Coulomb repulsion
https://scipost.org/SciPostPhys.19.3.067
doi: https://doi.org/10.21468/SciPostPhys.19.3.067



## Parameters

The parameter files in the `params` directory are used to control the system parameters.

### Temperature in K
`T <float>`
### Effective electron-electron interaction strength *g*
`phonon_coupling <float>`
### Debye frequency in meV
`omega_debye <float>`
### Fermi wavevector in sqrt(eV) - note that [hbar k_F / sqrt(m_e)] = sqrt(eV)
`k_F <float>`
### Scales the coulomb interaction, 1=physical reality, 0=turns off the Coulomb interaction
`coulomb_scaling <float>`
### Screening strength *lambda*
`screening <float>`
### Number of abscissae
`discretization_points <int>`
### Where to place the cutoff for the inner mesh (see the article)
`x_cut <float>`
### Name of the data dir
`output_folder <string>`


### Required externals
- Eigen https://eigen.tuxfamily.org/index.php?title=Main_Page or https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
- Boost https://www.boost.org/
- (optional) Intel MKL https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html - if present on the system, it will be used for Eigen. If not, the standard Eigen routines are employed.
- OpenMPI https://www.open-mpi.org/
- OpenMP https://www.openmp.org/
- CMake https://cmake.org/
- nlohmann/json.hpp https://github.com/nlohmann/json

### Required internals

- mrock/utility
- mrock/symbolic_operators

see https://github.com/majesticrock/PhdUtility


## Building

Build with `make`.
The Makefile handles the calls to cmake for you.
Without specification the project will be built for the local machine.
Specifying `cascadelake` or `icelake` will built for the specific CPU architecture.
These two are present on the compute cluster used to generate the data.


## Running the program

Before executing the program, make sure to build and run FermionCommute
https://github.com/majesticrock/FermionCommute/
such that you have the directory `../commutators/continuum/` filled with the results of the commutators.
`./exec.sh` will execute the program with the parameter file `params/params.config`.
For large scale computations, SLURM scripts are provided that use `params/cluster.config` and `params/ul_cluster.config`, respectively.
A few additional bash scripts are provided for ease of running multiple jobs.