# README #

### What is this repository for? ###

OptaDOS is a program for calculating core-electron and low-loss electron energy loss spectra (EELS) and optical spectra along with total-, projected- and joint-density of electronic states (DOS) from single-particle eigenenergies and dipole transition coefficients.

Energy-loss spectroscopy is an important tool for probing bonding within a material. Interpreting these spectra can be aided by first principles calculations. The spectra are generated from the eigenenergies through integration over the Brillouin zone. An important feature of this code is that this integration is performed using a choice of adaptive or linear extrapolation broadening methods which we have shown produces higher accuracy spectra than standard fixed-width Gaussian broadening. OptaDOS currently interfaces out-of-the-box with CASTEP and ONETEP and may be straightforwardly interfaced to any electronic structure code.

### How do I get set up? ###

It is reccomended to download the latest release of OptaDOS from the GitHub repository https://github.com/optados-developers/optados

Inside the top level optados/ directory are a number of subdirectories, `documents`, `examples`. The code may be compiled using the Makefile in the `optados` directory. 

The `SYSTEM`, `BUILD`, `COMMS_ARCH` and `PREFIX` flags must be set, either in the make.system, or from the command line (for example:
`make BUILD=fast`

No external libraries are required for serial execution.

*`SYSTEM:` Choose which compiler to use to make OptaDOS. The valid values are: g95 (default), gfortran, ifort, nag, pathscale, pgf90 and sun. Other compilers may be added manually by editing the make.system file.

*`BUILD:` Choose the level of optimisations required when making OptaDOS. The valid values are: fast (default): all optimisations or debug: no optimisations. All compiler warnings. Makes a binary containing full debugging information in operating system’s native format, including information suitable for analysing using the gprof profiler.

*`COMMS_ARCH:` Whether to compile for serial or parallel execution. The valid values are: serial (default) or mpi.

*`PREFIX:` Choose where to place the OptaDOS binary. The default is the OptaDOS directory.

### Who do I talk to? ###

* The Mailing List OPTADOS@JISCMAIL.AC.UK
