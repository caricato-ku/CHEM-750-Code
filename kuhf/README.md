 This is bootstrap restricted Hartree-Fock program written by the
 Spring 2015 class of Advanced Quantum Mechanics (CHEM 850) at the
 University of Kansas: Tal Aharon, Matthew Barclay, Allyssa Massie,
 Christopher Otolsky, Sijin Ren, Derek Rice, and Marco Caricato.

 The program is written in FORTRAN 08, it's self contained, and can
 be compiled with gfortran. The program is proudly inefficient, most
 likely contains bugs, and is intended for teaching use only. Anyone
 is welcome to use and modify it, but we do not technical support of
 any kind. Have fun with it!

 The integrals are computed with uncontracted s-type functions
 only. The exponents are read-in in a "alphas" file with the
 following format: The first row reports the number of elements, and
 the total number of functions. After each elements, the number of
 functions for that elements are reported. The program automatically
 discards the functions of elements that are not reported in the
 input geometry. Example:

 3, 16\
 H, 3\
        3.425250914\
 	      0.62391373\
        0.1688554040\
 C, 3\
       71.61683735\
       13.04509632\
       3.530512160\
 O, 10\
       0.1307093214D+03\
       0.2380886605D+02\
       0.6443608313D+01\
       0.5033151319D+01\
       0.9996722919D-01\
       0.1559162750D+00\
       0.1169596125D+01\
       0.3995128261D+00\
       0.6076837186D+00\
       0.7001154689D+00


 The geometry is reported in a "geometry" file in Cartesian
 coordinates in Angstrom. The first line reports the number of
 elements in the molecule. Example:

 3\
 O        0.000000    0.000000    0.110851\
 H        0.000000    0.783837   -0.443405\
 H        0.000000   -0.783837   -0.443405

