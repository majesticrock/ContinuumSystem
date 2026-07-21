# ContinuumSystem

Code used to compute the mean-field data and collective excitation spectra presented in 

Collective modes in superconductors including Coulomb repulsion
J. Althüser and G. S. Uhrig
https://doi.org/10.21468/SciPostPhys.19.3.067




### Requirements
- C++ 20 and a functioning compiler (tested with g++ 13.3.0 on WSL and g++ 11.5.0 on Red Hat)
- Eigen (tested with version 3.4.1) https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
- Boost (tested with version 1.78.0) https://www.boost.org/
- OpenMPI (tested with version 4.1.1) https://www.open-mpi.org/
- OpenMP https://www.openmp.org/
- nlohmann/json.hpp (tested with version 3.11.3) https://github.com/nlohmann/json
- [recommended] CMake 3.30 or newer (tested with version 3.31.8)  https://cmake.org/
- [optional] Intel MKL for BLAS and LAPACK (tested with version 2025.2.1) https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html 

- the mrock library located in `../../PhdUtility/`.
- Before executing the program, make sure to build and run `FermionCommute` in `../FermionCommute/`.



## Parameters

The parameter files in the `params` directory are used to control the system parameters.

| Option | Type / Values | Description |
|---|---:|---|
| `T` | `<float>` | Temperature in K. |
| `phonon_coupling` | `<float>` | Effective electron–electron interaction strength \(g\). |
| `omega_debye` | `<float>` | Debye frequency in meV. |
| `k_F` | `<float>` | Fermi wavevector in \(\sqrt{\mathrm{eV}}\). Note: \([\hbar k_F/\sqrt{m_e}] = \sqrt{\mathrm{eV}}\). |
| `coulomb_scaling` | `<float>` | Scales the Coulomb interaction (1 = physical, 0 = Coulomb off). |
| `screening` | `<float>` | Screening strength \(\lambda\). |
| `discretization_points` | `<int>` | Number of abscissae. |
| `x_cut` | `<float>` | Where to place the cutoff for the inner mesh (see the article). |
| `output_folder` | `<string>` | Name of the data directory. |



## Building

Build with `make`.
The Makefile handles the calls to cmake for you.
Without specification the project will be built for the local machine (-march=native).
There are additionally the targets `icelake` and `cascadelake`, which will built for the corresponding CPU architecture.


## Running the program

Before executing the program, make sure to build and run `FermionCommute` in `../FermionCommute/`.
such that you have the directory `../commutators/hubbard/` filled with the results of the commutators.
Then, (after building of course), you may create or edit a parameter file in the `params` directory.
Run the program with `./path/to/executable path/to/param/file.config`.
By default, the executable will be located in `./build/default/`.
The result will be saved into `../../data/continuum/<output_folder in the parameter file>/<string generated from the model parameters>/<filename (e.g. resolvents.json.gz)>`.



## Testing

`make test` will build and run the program with a set of parameters listed in the `tests` directory.
It will then call the plot script in the same dir to visualize the test data.
As before, it requires a previously completed run of `FermionCommute`, which must save the commutators to `../commutators/hubbard/`.
The program will save the simulation data and consequently plot it.

The plots will be saved to `build/default/<plot_name>.pdf`.
For doing so, python with matplotlib, numpy, and pandas is required.
Moreover, the python modules of `../../PhdUtility/python` are needed to evaluate the continued fractions.

### Expected results
#### No Coulomb
Phase peak at ω=0
Higgs peak at ω=2Δ

#### Normal screening
Phase peak at ω>0  (around 5.5 meV)
Higgs peak at ω<2Δ (around 7.5 meV)

#### Weak screening
No visible phase peak 
Higgs peak at ω=2Δ 

#### Strong attraction
Phase peak at ω=0
Multiple peaks at ω<2Δ