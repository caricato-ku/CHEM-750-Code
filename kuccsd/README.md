  This is bootstrap unrestricted CCSD program written by the Spring
  2015 class of Advanced Quantum Mechanics (CHEM 850) at the
  University of Kansas: Tal Aharon, Matthew Barclay, Sijin Ren, and
  Marco Caricato.

  This program is designed to accept previously calculated
  Hartree-Fock two-electron integrals and apply the coupled-cluster
  post-SCF method.

  This program will be carried out in a number of steps:\
    (1) Reading in the previously calculated integrals\
    (2) Calculating an initial guess of the excitation Amplitudes
          defined solely as the corresponding two-electron integrals\
    (3) Iteratively solving for the excitation Amplitudes and
    incorporating them into the final coupled-cluster energy

  This first block creates a Module for Global variables, allowing us
  to separate our code into a number of Subroutines that each
  calculate various Terms of the coupled-cluster Amplitude equation.

 Input file:\
 Method = choice of method (CCD or CCSD)\
 EpsI = occupied orbital energies\
 EpsA = virtual orbital energies\
 abcd = <ab||cd> integrals (in that index order)\
 iabc = <ia||bc> integrals (in that index order)\
 iajb = <ia||jb> integrals (in that index order)\
 ijab = <ij||ab> integrals (in that index order)\
 ijka = <ij||ka> integrals (in that index order)\
 ijkl = <ij||kl> integrals (in that index order)\

 The name of the files is generic. If no file name is specified, the
 default names are: EpsA, EpsI, Iabcd, Iiabc, Iiajb, Iijab, Iijka, Iijkl.

