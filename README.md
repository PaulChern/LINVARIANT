# LINVARIANT
**LINVARIANT** is a Universal model generator, in the field of computational condensed matter physics. <br />
**INVARIANT** are such objects though we construct with the help of coordinate systems turn out to be independent of coordinates. <br />
**L** is to memorize famous physicist Lev Davidovich Landau (22 January 1908 – 1 April 1968). <br />
- Utilizing the Group Theory, it mathematically modeling physics systems such as Lattice, Electron, Spin and their coupling systems.
- LINVARIANT generates symmetry adapted microscopic or phenomenological models.
- LINVARIANT creates DFT training sets for the model fitting.
- LINVARIANT has numerical solvers running on the model, such as MC, MD, Exact diagnalization, Minimization, and et al..
- For large scale calculations, LINVARIANT exports the symbolic models into FORTRAN modules.
## Features:
### models:
- Displacement and magnetic modes analysis
- INVARIANT generator (couplings among distortion, strain and magnetic)
- Modulate modes in supercell, generate initial Domain walls (DWs) structures for Density Functional Theory (DFT) codes
- Landau model builder, Generate DFT training sets by phonon and magnetic frozen-in.
- Heisenberg model with DM model builder
- Tight-Binding model builder
### solvers:
- (1) Finite Element Method (FEM), (2) Minimization, (3) molecular dynamics (MD), (4) Monte Carlo (MC), and (5) Finite Differences nonlinear solver on large scale continuous model
- Parallel tempering algorithm is available with both MC and MD
### fitting:
- Basis (ionic): phonon/irreducible representation/atomistic basis
- Basis (electronic): pesudo-atomic/Wannier basis
- structures searching by energy invariants
- machine learning of energy invariants and their energy potential surfaces
- supervised model fitting
### Auxiliary:
- Write Fortran using mathematica
- interface to VASP, Quantum Espresso, and OpenMX
- interface to WANNIER90
- mpi and openmp parallelization
- dynamics under external electric field
- Jij of Heisenberg model from DFT by Liechtenstein formalism
- Fij (force constants) from tight-binding models
- Electron/phonon bands unfolding
- phonon/magnon calculations from DFT input
- X ray diffraction simulation
- Nudged Elastic Bands (NEB) and Growing String Method (GSM) to explore the phase transition, dynamics, and domain wall structures
- Mollwide projection
### Examples:
- Boracite, Perovskite (To be added: Spinel, Rutile, Pyrochlore)
## Todo:
- implement DFTB
- implement LLG dynamics
- implement k dot p model builder
- add exact diagonalization solver
- adding Atomistic Green's Function (AGF) method to study scattering and ultrafast non-equilibrium dynamics
- including electron phonon coupling (EPC) by fitting phonon dependent Tight-Binding (TB) model from DFT molecular dynamics (MD).
## Authors
* **Peng Chen** - peng.chen.iphy@gmail.com
* **Hongjian Zhao** - solidstatezhao@gmail.com
* **Sergey Artyukhin** - sergey.artyukhin@iit.it
* **Laurent Bellaiche** - laurent@uark.edu   <br />
See also the list of [contributors](https://github.com/PaulChern/LINVARIANT/contributors) who participated in this project.
