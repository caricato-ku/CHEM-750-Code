! This is bootstrap restricted Hartree-Fock program written by the
! Spring 2015 class of Advanced Quantum Mechanics (CHEM 850) at the
! University of Kansas: Tal Aharon, Matthew Barclay, Allyssa Massie,
! Christopher Otolsky, Sijin Ren, Derek Rice, and Marco Caricato.
!
! The program is written in FORTRAN 08, it's self contained, and can
! be compiled with gfortran. The program is proudly inefficient, most
! likely contains bugs, and is intended for teaching use only. Anyone
! is welcome to use and modify it, but we do not technical support of
! any kind. Have fun with it!
!
! The integrals are computed with uncontracted s-type functions
! only. The exponents are read-in in a "alphas" file with the
! following format: The first row reports the number of elements, and
! the total number of functions. After each elements, the number of
! functions for that elements are reported. The program automatically
! discards the functions of elements that are not reported in the
! input geometry. Example:
!
! 3, 16
! H, 3
!        3.425250914
! 	 0.62391373
!        0.1688554040
! C, 3
!       71.61683735
!       13.04509632
!       3.530512160
! O, 10
!       0.1307093214D+03
!       0.2380886605D+02
!       0.6443608313D+01
!       0.5033151319D+01 
!       0.9996722919D-01
!       0.1559162750D+00
!       0.1169596125D+01
!       0.3995128261D+00
!       0.6076837186D+00
!       0.7001154689D+00
!
!
! The geometry is reported in a "geometry" file in Cartesian
! coordinates in Angstrom. The first line reports the number of
! elements in the molecule. Example:
!
! 3
! O        0.000000    0.000000    0.110851
! H        0.000000    0.783837   -0.443405
! H        0.000000   -0.783837   -0.443405
!
!
!**********************************************************************
! Integrals module. 
! Here it also reads the basis set file (alphas) and geometry file (geometry)
Module KUHF_Integrals
CONTAINS
Subroutine Integrals(r,Nele,Vnn,S,Core)
  !this is the master program for calculating integals
  !initalizing  a through z as integers
  implicit integer (a-z) 
  !
  !creation of named integers, amazingchris is a counter 
  integer atoms, amazingchris, Nele
  !
  !symbol holds the element from reading in atomic positions, elements
  !holds it in for gaussians
  character(3), allocatable, dimension(:) :: symbol, elements 
  !
  !xs, ys, and zs, are the atomic positions. ls, ly, lz are the
  !gaussian's positions. tal holds all the exponents for each element,
  !exponents holds the exponents based on the coordinates (with respect
  !to the lx, ly, lz array)
  double precision, allocatable, dimension(:) :: xs, ys, zs,lx,ly,lz,&
       tal, exponents
  !
  !allocation of arrays with two indicies for distances between
  !gaussians, S(overlap), kinetic, and one electron potential
  double precision, allocatable, dimension(:,:) :: distances, S, Core, oneKE, oneV
  ! 
  !initalizing some temporary variables                
  double precision angstrom, deltax, deltay, deltaz, pi, temp, hold, distRc, &
       distTC, rij, Vnn
  double precision, allocatable, dimension(:,:,:,:) :: twoV
  !
  !inital positions for the product of two gaussians, Rp is for one
  !electron potential, Rik and Rjl are for two electron potential
  double precision, dimension (1:3) :: Rp, Rik, Rjl
  !
  !array of integers. numgau holds how many gaussians are at each
  !position corresponding to the atomic positions, charge is the atomic
  !charge, pos tells us where in the tal array the gaussians of each
  !element are
  integer, allocatable, dimension(:) :: numgau, charge, pos
  !
  !conversion from angstrom to bohr radii
  angstrom = .52917721092D0
  !This is Pie
  pi = 4.*atan(1.) 
  !
  !The first line of alphas is the # of elements. Then the symbol on the
  !next line, followed by how many gaussians are associated with
  !it. Then the actual gaussians
  open(15,file='alphas',status='old')
  read(15,*) f, t
  allocate(elements(f))
  allocate(numgau(f))
  allocate(tal(t))
  !  
  !x will be a counter, because the loop resets and the array needs to
  !continue (it needs to go as the sum of all the j's) 
  !
  !put the exponents
  !that correspond to the atomic symbol (the alphas file holds each
  !element and the exponents of the s-type gaussians that go with those
  !elements) into the tal array
  x = 1
  do i = 1, f 
     read(15,*) elements(i), numgau(i)
     do j = 1, numgau(i)
        read(15,*) tal(x)
        x = x + 1
     end do
  end do
  close(15)
  ! print*,'elements',elements
  ! print*,'numgau',numgau
  ! print*,'functions',tal
  !
  !open file with atom symbol and postions. assumes the first line of
  !the file is the number of atoms in the file (and consequently the
  !number of lines remaining)
  open(11,file='geometry',status='old')
  read(11,*) atoms
  !allocate arrays to hold what is read from the file with the atomic positions
  allocate (symbol(atoms))
  allocate (xs(atoms))
  allocate (ys(atoms))
  allocate (zs(atoms))
  !
  !allocate the array which will hold the atomic charges (for
  !potential integral)
  allocate (charge(atoms))
  !	
  !read in the atomic symbol and x, y, and z coordinates
  do i = 1, atoms 
     read(11,*)symbol(i),xs(i),ys(i),zs(i)
  end do
  close(11)
  !
  !MC use geometry in bohr for now
  !converts coordinates from angstroms to bohr radii
  do i =1, atoms
     xs(i) = xs(i) / angstrom
     ys(i) = ys(i) / angstrom
     zs(i) = zs(i) / angstrom
  end do
  !
  !gets the atomic charge of each atom and puts it in an array, both for
  !figuring out which gaussians to use and for the potential integral
  Nele = 0
  do i = 1, atoms
     if (symbol(i) == "H") then
        charge(i) = 1
        Nele = Nele + 1
     else if (symbol (i) == "He") then
        charge(i) = 2
        Nele = Nele + 2
     else if (symbol (i) == "Li") then
        charge(i) = 3
        Nele = Nele + 3
     else if (symbol (i) == "Be") then
        charge(i) = 4
        Nele = Nele + 4
     else if (symbol (i) == "B") then
        charge(i) = 5
        Nele = Nele + 5
     else if (symbol (i) == "C") then
        charge(i) = 6
        Nele = Nele + 6
     else if (symbol (i) == "N") then
        charge(i) = 7
        Nele = Nele + 7
     else if (symbol (i) == "O") then
        charge(i) = 8
        Nele = Nele + 8
     else if (symbol (i) == "F") then
        charge(i) = 9
        Nele = Nele + 9
     else if (symbol (i) == "Ne") then
        charge(i) = 10
        Nele = Nele + 10
     end if
  end do
  !
  r = 0
  !how long will the lx, ly, and lz array be, and more importantly,
  !how many total gaussians will be used for the molecule
  do i = 1, atoms
     do j = 1, f
        !numgau stores the number of gaussians that go with each
        !element. when the loop finds an atomic symbol in the symbols
        !array, it adds the number of gaussians that go with that
        !element
        if (symbol(i) == elements(j)) then
           r = r + numgau(j)
        end if
     end do
  end do
  print*,"N basis=",r
  print*,"N electrons=",Nele
  !
  !allocate some arrays and matrices, broskoi
  allocate (S(r, r))
  allocate (oneKE(r, r))
  allocate (oneV(r, r))
  allocate (twoV(r, r, r, r))
  allocate (lx(r))
  allocate (ly(r))
  allocate (lz(r))
  allocate (distances(r, r))
  allocate (exponents(r))
  allocate (pos(f))
  !
  !position tells us where in the tal array we can find the exponents
  !that should go with the gaussians of individual elements
  pos(1) = 0
  do g = 2, f
     pos(g) = pos(g-1) + numgau(g-1)
  end do
  !
  !this makes array with the position of each gaussian (lx, ly, lz),
  !and then puts an the exponent that goes with that position in the
  !exponents array
  amazingchris = 1
  do i = 1, atoms
     do k = 1, f
        if (symbol(i) == elements(k)) then
           do j = 1, numgau(k)
              lx(amazingchris) = xs(i)
              ly(amazingchris) = ys(i)
              lz(amazingchris) = zs(i)
              !pos tells us where to find our element, and pos(k)
              !corresponds to a specific element. we add j because the
              !same element can have multiple different gaussians with
              !different exponents
              exponents(amazingchris) = tal(pos(k)+j)
              !chris likes to count
              amazingchris = amazingchris + 1
           end do
        end if
     end do
  end do
  !MC
  ! Compute nuc-nuc repulsion
  Vnn = 0.d0
  Do i = 1, atoms-1
     Do j = i+1, atoms
        rij = (xs(i)-xs(j))**2 + (ys(i)-ys(j))**2 + (zs(i)-zs(j))**2  
        Vnn = Vnn + charge(i)*charge(j)/sqrt(rij)
     end Do
  end Do
  !
  !fill distances with 0's for easier visualization
  distances = 0.0
  ! do i = 1, r
  !    do j = 1, r
  !       distances(i,j) = 0
  !    end do
  ! end do
  !
  !filling more matrices with 0's, just to initialize
  S = 0
  oneKE = 0
  oneV = 0
  ! do i = 1, r 
  !    do j = 1, r
  !       S(i,j) = 0
  !       oneKE(i,j) = 0
  !       oneV(i,j) = 0
  !    end do
  ! end do
  !
  !calculate the distance between gaussians
  do i = 1, r
     do j = 1, r
        distances(i,j) = sqrt((lx(i)-lx(j))**2 &
             +(ly(i)-ly(j))**2 &
             +(lz(i)-lz(j))**2)
     enddo
  enddo
  ! print*,'exponents'
  ! do i =1, r
  !    print*,i,lx(i),ly(i),lz(i),exponents(i)
  ! end do
  !
  do i = 1, r 
     do j = 1, r
        !calculate the overlap and one electron kinetic energy integrals
        call overlap(exponents(i),exponents(j),distances(i,j),S(i,j),&
             &oneKE(i,j))         
        do k = 1, atoms
           !get the position of the "Rp" gaussian, which is the one
           !that is the product of the two gaussian for the one
           !electron nuclear attraction integral
           !
           !get the x cooridnate of the Rp gaussian (product of
           !gaussians i and j)
           Rp(1) = (exponents(i) * lx(i) + exponents(j) * lx(j))/& 
                &(exponents(i)+exponents(j))

           !get the y cooridnate of the Rp gaussian (product of
           !gaussians i and j)
           Rp(2) = (exponents(i) * ly(i) + exponents(j) * ly(j))/& 
                &(exponents(i)+exponents(j))
           !get the z cooridnate of the Rp gaussian (product of
           !gaussians i and j)
           Rp(3) = (exponents(i) * lz(i) + exponents(j) * lz(j))/& 
                &(exponents(i)+exponents(j))
           !calculates the distance between the Rp gaussian and the
           !charge (in this case the atom)
           distRc = sqrt((Rp(1) - xs(k))**2 + (Rp(2) - ys(k))**2 +&
                & (Rp(3) - zs(k))**2) 
           !calculate the one electron potential
           call onepot(exponents(i), exponents(j), distances(i,j)&
                &, distRc, charge(k), oneV(i,j))
        end do
     end do
  end do
  do l = 1, r
     do k = 1, r
        do j = 1, r
           do i = 1, r
              !get the x cooridnate of the Rik gaussian (product of
              !gaussians i and k)
              Rik(1) = (exponents(i) * lx(i) + exponents(k) * lx(k))/& 
                   &(exponents(i)+exponents(k))
              !get the y cooridnate of the Rik gaussian (product of
              !gaussians i and k)
              Rik(2) = (exponents(i) * ly(i) + exponents(k) * ly(k))/& 
                   &(exponents(i)+exponents(k))
              !get the z cooridnate of the Rik gaussian (product of
              !gaussians i and k)
              Rik(3) = (exponents(i) * lz(i) + exponents(k) * lz(k))/& 
                   &(exponents(i)+exponents(k))
              !get the x cooridnate of the Rjl gaussian (product of
              !gaussians j and l)
              Rjl(1) = (exponents(j) * lx(j) + exponents(l) * lx(l))/&
                   &(exponents(j)+exponents(l))
              !get the y cooridnate of the Rjl gaussian (product of
              !gaussians j and l)
              Rjl(2) = (exponents(j) * ly(j) + exponents(l) * ly(l))/& 
                   &(exponents(j)+exponents(l))
              !get the z cooridnate of the Rjl gaussian (product of
              !gaussians j and l)
              Rjl(3) = (exponents(j) * lz(j) + exponents(l) * lz(l))/& 
                   &(exponents(j)+exponents(l))
              !we get the distances between the two gaussians that are
              !the products of the four gaussians
              distTC = sqrt((Rik(1) - Rjl(1))**2 + (Rik(2) - Rjl(2))**2&
                   & + (Rik(3) - Rjl(3))**2)
              !calculates the two electron integrals
              call electron(exponents(i),exponents(j),exponents(k),&
                   &exponents(l),distances(i,k),distances(j,l),distTC,&
                   &twoV(i,j,k,l))
           end do
        end do
     end do
  end do
  !  
  ! writes all of the output matrices to seperate text files, for
  ! communication with the other Subroutines for the master program
  ! open(unit=10,file='s.txt')
  ! write(10,*) r
  ! write(10,*) S
  ! close(10)
  ! open(unit=12,file='kinetic.txt')
  ! write(12,*) r
  ! write(12,*) oneKE
  ! close(12)
  ! open(unit=13,file='oneelectronpotential.txt')
  ! write(13,*) r
  ! write(13,*) oneV
  ! close(13)
  open(unit=14,file='twoelectronpotential.txt')
  write(14,*) r
  write(14,*) twoV
  close(14)
  ! do l = 1, r
  !    do k = 1, r
  !       do j = 1, r
  !          do i = 1, r
  !             print*, l,k,j,i, twoV(i,k,j,l)
  !          end do
  !       end do
  !    end do
  ! end do
  ! print*, "2,2,2,2", twoV(2,2,2,2)
  ! print*, "2,2,2,1", twoV(2,2,2,1)
  ! print*, "2,2,1,1", twoV(2,2,1,1)
  ! print*, "2,1,2,1", twoV(2,1,2,1)
  !
  ! Do i = 1, r
  !    Do j = 1, r
  !       Core(i,j) = oneKE(i,j) + oneV(i,j)
  !       S(i,j) = S(i,j)
  !    end Do
  ! end Do
  ! open(unit=15,file='CORE.txt')
  ! write(15,*) r
  ! write(15,*) Core
  ! close(15)
  ! print*, "Kinetic"
  ! Do i = 1, r
  !    print*, (oneKE(i,j),j=1,r)
  ! end Do
  ! print*, "Potential"
  ! Do i = 1, r
  !    print*, (oneV(i,j),j=1,r)
  ! end Do
  Core = oneKE + oneV
  Deallocate(elements,numgau,tal,symbol,xs,ys,zs,charge,oneKE,oneV,&
       &twoV,lx,ly,lz,distances,exponents,pos)
  !  Deallocate(smat)
  !
End Subroutine Integrals
!
!***************************************************************
Subroutine fosho(arg, retval)
  implicit double precision (A-H, O-Z)
  double precision, intent(in) :: arg
  double precision, intent(out) :: retval
  double precision pi, morepi
  pi = 4.*atan(1.) !This is pie
  !assymtotic condition, as x => 0, erf(x) => 1
  if (abs(arg) .lt. 1.0D-6) then
     retval = 1.0D0 - arg
     !calls the error function
  else 
     call derf(sqrt(arg), morepi)
     !returns the value in the Fo function
     retval = 0.5D0 * sqrt(pi/arg) * morepi 
  end if
  !
end Subroutine fosho
!
!***************************************************************
Subroutine derf(arg, retval) 
  !polynomial approximation of the error function
  ! This is taken from "Modern Quantum Chemistry" By Szabo and Ostlund
  implicit double precision (A-H, O-Z)
  double precision, intent(in) :: arg
  double precision, intent(out) :: retval
  Dimension A(5)
  double precision TN, poly
  data P/0.3275911D0/
  Data A/0.254829592D0, -0.284496736D0, 1.421413741D0, -1.453152027D0,&
       & 1.061405429D0/
  T = 1.0D0/(1.0D0 + P*arg)
  TN = T
  poly = A(1) * TN
  do I = 2, 5
     TN = TN * T
     poly = poly + A(i) * TN
  end do
  retval = 1.0D0 - poly * exp(-arg*arg)
end Subroutine derf
!
!***************************************************************
Subroutine overlap (alpha1, alpha2, dist, retval,retval2)
  implicit none
  double precision, intent(in) ::  alpha1, alpha2, dist
  double precision, intent(out) :: retval,retval2
  double precision :: pi, ss, ke
  !This is pie
  pi = 4.*atan(1.) 
  !calculates the overlap
  !calculation of s-s type overlap
  ss = (((8.*alpha1**3)/(pi**3))**.25)&
       &*(((8.*alpha2**3)/(pi**3))**.25)&
       &*((pi/(alpha1+alpha2))**1.5)&
       &*exp(((-alpha1*alpha2)/(alpha1+alpha2))*(dist**2))
  !calculates kinetic energy
  !calculation of the S-S type kinetic energy
  ke = (((8.*alpha1**3)/(pi**3))**.25)&
       &*(((8.*alpha2**3)/(pi**3))**.25)& 
       &*(alpha1*alpha2)/(alpha1+alpha2)&
       &*(3-2*(((alpha1*alpha2)/(alpha1+alpha2))*(dist**2)))*&
       &((pi/(alpha1+alpha2))**1.5)&
       &*exp(((-alpha1*alpha2)/(alpha1+alpha2))*(dist**2))
  !send them back
  retval = ss 
  retval2 = ke
end Subroutine overlap
!
!***************************************************************
Subroutine onepot(alpha1, alpha2, distAB, distCP, charge, retval)
  implicit double precision (a-e)
  integer, intent(in) :: charge
  double precision, intent(in) ::  alpha1, alpha2, distAB, distCP
  double precision, intent(out) :: retval
  double precision :: pi, ss, temp, mid
  !This is pie
  pi = 4.*atan(1.) 
  !temp is what we pass into the error function
  temp = (alpha1 + alpha2)*(distCP**2)
  !call Subroutine to calculate the error functon
  call fosho(temp, mid)
  !the actual calculation of the one electron potential integral
  ss = (((8.*alpha1**3)/(pi**3))**.25)&
       &*(((8.*alpha2**3)/(pi**3))**.25)*((-2.0*pi)&
       &/(alpha1 + alpha2))*charge&
       &*exp(((-alpha1*alpha2)/(alpha1+alpha2))*(distAB**2)) * mid
  !this is a sum because the values for each different charge but same
  !gaussians are summed
  retval = retval + ss
end Subroutine onepot
!
!***************************************************************
Subroutine electron(a1, a2, a3, a4, distAB, distCD, distTC, retval)
  implicit none
  double precision, intent(in) ::  a1, a2, a3, a4, distAB, distCD, distTC
  double precision, intent(out) :: retval
  double precision :: pi, pi3, pi5, ss, temp, mid, norm, fact, ee
  !This is pie
  pi = 4.*atan(1.) 
  pi3 = pi*pi*pi
  pi5 = pi*pi*pi*pi*pi
  !temp is what we actually pass into the error function
  temp = (a1+a3)*(a2+a4)/(a1+a2+a3+a4)*(distTC**2)
  !calculate the error function
  call fosho(temp, mid)
  !calculate the two electron integral
  !MC improve readability
  norm = a1*a2*a3*a4
  norm = 8.d0*sqrt(sqrt(norm**3))/pi3
  fact = (a1+a3)*(a2+a4)*sqrt(a1+a2+a3+a4)
  fact = 2.d0*sqrt(pi5)/fact
  ee = a1*a3*(distAB**2)/(a1+a3) + a2*a4*(distCD**2)/(a2+a4)
  ss = norm*fact*exp(-ee)*mid
  ! ss =((((8.*a1**3)/(pi**3))**.25)*(((8.*a2**3)/(pi**3))**.25)&
  !      &*(((8.*a3**3)/(pi**3))**.25)*(((8.*a4**3)/(pi**3))**.25))&
  !      &*2*(pi**2.5)/((a1+a3)*(a2+a4)*((a1+a2+a3+a4)**.5))&
  !      &*exp((((-a1*a3)/(a1+a3))*(distAB**2))-(((a2*a4)/(a2+a4))*(distCD**2)))&
  !      &*mid
  retval = ss
  !  print*, distAB,distCD,distTC,exp(-ee),temp,mid,ss
end Subroutine electron

End Module KUHF_Integrals

!
!**********************************************************************
Subroutine SDiag(N,S,U)
  !Program will diagonalize overlap matrix S
  implicit none
  integer :: i, j, k,n,o,p,count,maxcycle
!  double precision, allocatable, dimension (:,:) :: S, U        !S is overlap U
  double precision, dimension (N,N) :: S, U
  double precision ode, oav, beta, c, si, t, sic, csi,pivot         !Defined later
  double precision, parameter :: thresh = 1.d-12
  ! open (13, file='s.txt')         !Opens the overlap matrix calc by Chris and Tal
  ! read (13,*) n                 !First line should give matrix size n
  ! allocate (S(n,n))             !Builds S array with n dimensions from file
  ! read (13,*) S                 !Fills in matrix S
  ! allocate (U(n,n))
  ! print*,'SDiag0'
  ! do i = 1, n
  !    print*, (S(i,j),j=1,n)
  ! end do
  U = 0.0                        !Initializes the empty identity matrix
  do i=1,n                       !After rotations, this will become unity matrix
     U(i,i)=1.0              !Fills in 1s down the diagonal
  end do
  !calculates squared sum of off diagonal elements variable ode
  !The program will terminate when these are below a certain threshold that will
  !be defined laer- once they;re less thatn 1.0e-10
  ode = 0.0
  do j=1,n
     do i=1,n
        if (i.ne.j) ode = ode + S(i,j)**2 !ignores diagonal elements
     end do
  end do
  !Calculates average of off-diagonal elements divided by 2 variable is oav
  oav = 0.5*ode/dble(n*n)
  !MC
  !determine first pivot here. This avoids awkward case where S is already diagonal
  pivot = 0.0
  do j= 1,n-1
     do i=j+1,n
        if (abs(S(i,j)).gt.pivot) then
           o = i
           p = j
           pivot = abs(S(i,j))
        end If
     end do
  end do
  !  pivot = .00001                  !Starting the pivot value before loop
  count = 1
  maxcycle = 1000
!  print*,'pivot-0', pivot,o,p,n
  do while (pivot >= thresh .and. count <= maxcycle )      !goes while sum of ode is larger than
     ! pivot = 0.0               !Starts the cycle
     ! do j= 1,n
     !    do i=1,n
     !       if (i.eq.j) cycle               !skips diagonal
     !       !this finds largest off diagonal elements and rotates it to get rid of it
     !       !MC
     !       if (abs(S(i,j)).gt.pivot) then      
     !          !if (S(i,j).gt.pivot) then      
     !          o = i                         !This process is
     !          p = j                         !to get rid of
     !          !MC
     !          pivot = abs(S(i,j))
     !          ! pivot = S(i,j)
     !       end if                !off diag elems
     !       !end if
     !    end do
     ! end do
     i = o
     j = p
     ode = ode- 2.0*S(j,i)**2
     oav = 0.5*ode/dble(n*n)
     !Next step calcs elems of Givens Matrix- beta,t,c,si-parameters for rotation
     beta = ((S(j,j)-S(i,i))/(2.0*S(j,i)))
     t = 0.5 * beta/SQRT(1.0 + beta**2)
     c = sqrt(max(0.5-t,0.0))        !cos
     si = sqrt(max(0.5+t,0.0))       !sin
     !Next step recalculates i and j
     do k=1,n
        csi = c*S(i,k)+si*S(j,k)  !cos*sin
        sic = -si*S(i,k)+c*S(j,k) !sin*cos
        S(i,k) = csi
        S(j,k) = sic
        !print*,k,i,j,csi,sic
     end do
     !Performs rotation to obtain transformed matrix S which will eventually be
     !diagonalized after various iterations
     !Matrix S is eigenvalues and U is eigenvectors-unitary matrix
     !These come from trig identities
     do k = 1,n
        csi = c*S(k,i)+si*S(k,j)
        sic = -si*S(k,i)+c*S(k,j)
        S(k,i) = csi
        S(k,j) = sic
        csi = c*U(k,i)+si*U(k,j)
        sic = -si*U(k,i)+c*U(k,j)
        U(k,i) = csi
        U(k,j) = sic
     end do
     pivot = 0.0               !Starts the cycle
     do j= 1,n
        do i=1,n
           if (i.eq.j) cycle               !skips diagonal
           !this finds largest off diagonal elements and rotates it to get rid of it
           !MC
           if (abs(S(i,j)).gt.pivot) then      
              !if (S(i,j).gt.pivot) then      
              o = i                         !This process is
              p = j                         !to get rid of
              !MC
              pivot = abs(S(i,j))
              ! pivot = S(i,j)
           end if                !off diag elems
           !end if
        end do
     end do
!     print *, count,pivot
     count = count + 1
  end do
  !end do
  !
  ! open(14, file='DiagS.txt')      !Creates new file
  ! write(14,*) n                 !writes dimensions of Matrix
  ! write(14,*) S                 !Writes out diagonalized s matrix
  ! 
  ! open(22, file='Unitary.txt')    !Creates new file for unitary matrix
  ! write(22,*) n                 !Writes dimensions
  ! write(22,*) U                 !writes out unitary matrix
  ! do i=1,n-1
  !    do j=i+1,n
  !       print *, i,j, 'Sij', S(i,j)
  !       print*, "test"
  !    end do
  ! end do
  ! 
  ! do i=1,n
  !    print *, i,'Sii', S(i,i)
  !    print*, "test2"
  ! end do
  !MC
  ! print*,'SDiag'
  ! do i = 1, n
  !    print*, (S(i,j),j=1,n)
  ! end do
End Subroutine SDiag
!
!**********************************************************************
Subroutine TransformationMatrix(N,S,X)   
  !Program will create transformation Matrix X
  ! for canonical orthogonalization
  implicit none
  integer :: i, j, n
  double precision, dimension (N,N) :: S, X 
!  double precision, allocatable, dimension (:,:) :: S, X, U 
  !
  ! open (12, file='DiagS.txt')       !Opens the diagonal overlap matrix d
  ! read (12,*) n             !First line should give matrix size n
  ! allocate (S(n,n))         !Builds s array with n dimensions
  ! read (12,*) S             !Fills in matrix S that is diagonalized
  ! !
  ! open (22, file='Unitary.txt')   
  ! read (22,*) n                 !First line should give matrix size n
  ! allocate (U(n,n))             !Builds s array with n dimensions
  ! read (22,*) U                 !Fills in unitary matrix
  ! allocate (X(n,n))                       !allocates tranform matrix
  ! FormX 
  ! eq 3.169 of Szabo
  do j=1,n
     do i=1,n
        !MC
        ! if (S(i,j)== 0) cycle            !To avoid diviing by zero
        ! X(i,j) = U(i,j)*(1/(SQRT(S(i,j))))      !Comes from eq 3.169
        X(i,j) = X(i,j)/sqrt(S(j,j))
     end do
  end do
  ! open(32, file='XMatrix.txt') !Creates new file
  ! write(32,*) n                 !writes dimensions of Matrix
  ! write(32,*) X                 !Writes out transform matrix X
  !
End Subroutine TransformationMatrix
!
!**********************************************************************
Subroutine Gmatrix(dim,G,P)
  !Program will calculate G matrix
  implicit none
  integer :: l, s, m, n, dim
  double precision, allocatable, dimension (:,:,:,:) :: I
  double precision, dimension (dim,dim) :: P, G
!  double precision, allocatable, dimension (:,:) :: P, G
!  double precision :: a, sum
  !
  open (1, file='twoelectronpotential.txt')                                 !open the 2e integral file
  read (1,*) dim                                                            !read in dimension of integral
  allocate (I(dim,dim,dim,dim))                                              !assign dimension for the 2e integral matrix
!  allocate (G(dim,dim))                                                     !assign dimension for the G matrix           
  read (1,*) I                                                              !read in and form 2e integral matrix
  close(1)
  G = 0.d0
  !
  ! open (2, file='P.txt')                                                    !open the P matrix file
  ! read (2,*) dim                                                            !read in dimension of integral
  ! allocate (P(dim,dim))                                                     !assign dimension for the P matrix
  ! read (2,*) P                                                              !read in and form P matrix
  !do loops that calculate the elements of G matrix
  do n=1,dim                                                                
     do m=1,dim
        do l=1,dim
           do s=1,dim
              G(m,n) = G(m,n) + P(l,s)*(I(m,s,n,l)-0.5*I(m,s,l,n))
              ! a = P(l,s)*(I(m,s,n,l)-0.5*I(m,s,l,n))
              ! sum = sum + a
           end do
        end do
        ! G(m,n) = sum
        ! sum = 0
     end do
  end do
  Deallocate (I)
  !assign dimension for the 2e integral matrix
  ! !creat and open a new file for G the matrix
  ! open (3, file='Gmatrix.txt')                                             
  ! !write the dimention of the matrix in the first line of the G matrix file
  ! write (3,*) dim                                                          
  ! !write the matrix 
  ! write (3,*) G                                                            
  !
End Subroutine Gmatrix
! !
! !**********************************************************************
! Subroutine FockMatrix(NB,Core,F)
!   !Program will calculate Fock Matrix F
!   implicit none
!   integer :: i, j, n
! !  double precision, allocatable, dimension (:,:) :: T, V, G, F
!   double precision, allocatable, dimension (:,:) :: G
!   double precision, dimension (:,:) :: Core, F
!   ! open (11, file='Kinetic.txt')	!Opens the matrix of KE integrals 
!   ! read (11,*) n			!First line should give matrix size n
!   ! allocate (T(n,n))		!Builds T array with n dimensions from file
!   ! read (11,*) T			!Fills in matrix T
!   ! 
!   ! open (12, file='oneelectronpotential.txt')
!   ! read (12,*) n
!   ! allocate (V(n,n))
!   ! read (12,*) V                
!   open (51, file='Gmatrix.txt')	!Opens matrix G that Sijin calcs
!   read (51,*) n			!First line should give matrix dimensions
!   allocate (G(n,n))		!Builds G array with matrix size n
!   read (51,*) G			!Fills in matrix G
!   !  allocate (F(n,n))             !Makes F n by n matriz
!   do j=1,n
!      do i=1,n
!         !!Next line comes from eq 3.153 and 3.154 from Szabo book
!         !        F(i,j) = V(i,j) + T(i,j) + G(i,j)
!         F(i,j) = Core(i,j) + G(i,j)
!      end do			!This and next line will terminate do loops
!   end do
!   ! open(61, file='F.txt')	!Creates new file
!   ! write(61,*) n			!writes dimensions of Fock Matrix
!   ! write(61,*) F			!Writes out Fock Matrix F
!   !
! End Subroutine FockMatrix
!
!**********************************************************************
Subroutine transformedF(dim,F,X)
  !the program calculate tranfromed fock matrix
  implicit none                                        
  !  double precision, allocatable, dimension (:,:) :: F, X, Xtrans, transformF, M
  double precision, allocatable, dimension (:,:) :: Xtrans, M
  double precision, dimension (dim,dim) :: F, X
  integer :: dim
  allocate (M(dim,dim))                  
  allocate (Xtrans(dim,dim))             
  !calculate the transposed X matrix
  XTrans = 0.0
  M = 0.0
  Xtrans = Transpose(X)                  
  M = matmul(Xtrans, F)                  
  !calculate the transfored Fock matrix
  F = matmul(M, X)              
  Deallocate(XTrans,M)
  !
  ! open (1, file='F.txt')                 !open the file of fock matrix generated by other group      
  ! read (1,*) dim                         !read in the dimension of the fock matrix
  ! allocate (F(dim,dim))                  !assign dimension for the fock matrix
  ! allocate (transformF(dim,dim))         !assign dimension for the transform fock matrix
  ! allocate (M(dim,dim))                  !assign dimension for M matrix
  ! read (1,*) F                           !read in and form the fock matrix
  ! close(1)
  ! !
  ! open (2, file='X.txt')                 !open the file of X matrix generated by other group
  ! read (2,*) dim                         !read in the dimension of the X matrix
  ! allocate (X(dim,dim))                  !assign dimension for the X matrix
  ! allocate (Xtrans(dim,dim))             !assign dimension for the transposed X matrix
  ! read (2,*) X                           !read in and form X matrix
  ! Xtrans = TRANSPOSE(X)                  !calculate the transposed X matrix
  ! M = matmul(Xtrans, F)                  
  ! transformF = matmul(M, F)              !calculate the transfored Fock matrix
  ! !
  ! open (3, file='F1.txt')                 !open a new file for storing the transformed fock matrix
  ! write (3,*) dim                        !write in the dimension for the transformed fock matrix
  ! write (3,*) transformF                 !write in the fock matrix
  ! close(3)
  !
end Subroutine transformedF
!
!**********************************************************************
Subroutine FDiag(N,F,C)
  !This program will diagonalize the F matrix to give the eigenvalues, E and the
  !eigenvectors C. This uses the Jacobi method for diagonalization where a
  !rotation matrix is used to kill one off-diagonal term at a time until there are
  !no (or very small) off-diagonal terms left. The rotation matrix is continually
  !updated until one giant rotation matrix is built that diagonalizes F. This
  !rotation matrix is the C matrix.
  implicit none
  integer :: i, j, k, l, n, count, maxcycle
  double precision, dimension(N,N) :: F,C
!  double precision, allocatable :: F(:,:),C(:,:)
  double precision :: offsum, offavg
  double precision, parameter :: thresh = 1.d-20
  double precision :: ct, cs, co, s, temp1, temp2, t,y,step
  !
  !  open(20,file='F1.txt')
  ! thresh = 1.d-12
  ! read(20,*) n
  ! allocate(F(n,n))
  ! allocate(C(n,n))
  ! read(20,*) F
  ! close(20)
  !
  C = 0.0
  do i=1,n
     C(i,i) = 1.0
  end do
  !The C matrix starts as the identity matrix before any rotations are made.

  do i = 1,n-1
     do j=i+1,n
        if ((F(i,j)-F(j,i))>1.d-6) then
!           print *, 'Warning, Fock matrix is not symmetric. Forcing symmetry'
           F(i,j) = F(j,i)
        end if
     end do
  end do
  !The above can be used to force symmetry if needed.
  offsum = 0.0
  do i=1,n-1
     do j=i+1,n
        offsum = offsum + F(j,i)**2
     end do
  end do
  !Sums up the off diagonal terms
  offavg = 0.5*offsum/dble(n*n-n)
  !Calculates an average of the off diagonal so that later the program can focus
  !on larger terms
  count = 1
  maxcycle = 50
  !MC
  ! Add max cycles check
  do while (offsum.gt.thresh.and.count<=maxcycle)
     do i=1,n-1
        do j=i+1,n
           !This focuses on the larger terms by killing things above the average. This cuts
           !the number of steps down greatly.
           if (F(j,i)**2 <= offavg) cycle  
           !MC
           ! re-evaluate off-diagonal sum
           offsum = 0.0
           do k=1,n-1
              do l=k+1,n
                 offsum = offsum + F(l,k)**2
              end do
           end do
           !offsum = offsum - F(j,i)**2
           offavg = 0.5*offsum/dble(n*n-n)
           !Recalculate the sum and average here
           ct = (F(j,j)-F(i,i))/(2.0*F(j,i))
           cs = 0.5*ct/sqrt(1.0+ct**2)
           s = sqrt(max(0.5+cs,0.0))
           co = sqrt(max(0.5-cs,0.0))
           !As stated, this program will be performing rotations about some angle, theta,
           !that will specifically kill the chosen terms one at a time. ct is the
           !cot(2theta) which can be converted to sin and cos using trig identities. The
           !sin(s) and cos(co) calculated are the ones used in the rotation matrix. Theta
           !needs not be explicitly calculated since all that needs to be known is that the
           !rotation must be set to F(i,j)=0. The max keeps the rotations consistant.
           do k=1,n
              temp1 =  co*F(i,k)+s*F(j,k)
              temp2 = -s*F(i,k)+co*F(j,k)
              F(i,k) = temp1
              F(j,k) = temp2
           end do
           !Start by rotating the rows
           do k=1,n
              temp1 =  co*F(k,i)+s*F(k,j)
              temp2 = -s*F(k,i)+co*F(k,j)
              F(k,i) = temp1
              F(k,j) = temp2
              !Finish building th enew F matrix, which is slowly becoming the E diagonal
              !matrix
              temp1 =  co*C(k,i)+s*C(k,j)
              temp2 = -s*C(k,i)+co*C(k,j)
              C(k,i) = temp1
              C(k,j) = temp2
              !Build the C matrix
           end do
           step = step + 1
        end do
     end do
     count = count + 1
  end do
  ! open(23,file='Energy.txt')
  ! open(24,file='C1.txt')
  ! write(23,*) n
  ! write(24,*) n
  ! write(23,*) F
  ! write(24,*) C
  ! close(23)
  ! close(24)
  ! !MC
  ! do i = 1, n
  !    print*, (F(i,j),j=1,n)
  ! end do
  !
End Subroutine FDiag
!**********************************************************************
Subroutine Cmatrix(dim,X,C1)
  !the program calculate C matrix
  implicit none
  double precision, allocatable, dimension (:,:) :: C2
  double precision, dimension (dim,dim) :: X, C1
!  double precision, allocatable, dimension (:,:) :: X, C1, C2
  integer :: dim
  !
  allocate (C2(dim,dim))                     !assign dimension fo the C matrix
  ! open (1, file='C1.txt')                    !open C' matrix file
  ! read (1,*) dim                             !read in the dimension of C' matrix
  ! allocate (C1(dim,dim))                     !assign dimension for the C' matrix
  ! allocate (C2(dim,dim))                     !assign dimension fo the C matrix
  ! read (1,*) C1                              !read in and form C'matrix
  ! open (2, file='X.txt')                     !open X matrix file
  ! read (2,*) dim                             !read in the dimension of X matrix
  ! allocate (X(dim,dim))                      !assign dimension for the X matrix
  ! read (2,*) X                               !read in and form X matrix
  !
  C2 = matmul(X,C1)                          !calculate the C matrix
  C1 = C2
  Deallocate(C2)
  !
  ! open (3, file='C.txt')                     !open another file for storing C matrix
  ! write (3,*) dim                            !write in the dimension of C matrix
  ! write (3,*) C2                             !write in the C matrix
  !
End Subroutine Cmatrix
!
!**********************************************************************
Subroutine NewDensity(N,Nele,C,P)
  !Forms the new density matrix from C
  implicit none
  integer :: n,i,j,k,Nele
!  double precision, allocatable, dimension(:,:) :: C
  double precision, dimension(N,N) :: P,C
!   open(21,file='Cmatrix.txt')
!   !Openening the freshly calculated C matrix.
!   read(21,*) n
!   !First line of file should give n dimension size of array
!   allocate(C(n,n))
! !  allocate(P(n,n))
!   !Builds the size of the C matrix and P, the density, matrix
!   read(21,*) C
!   !Bring the previously calculated C matrix into the program
  ! print*, "Guess C",n
  ! Do i = 1, N
  !    print*, (C(i,j),j=1,n)
  ! end Do
  P = 0.d0
!  i = n/2
!  P = 2.0 * MatMul(C(:,1:i),Transpose(C(:,1:i)))
  do j=1,n
     do i=1,n
        do k=1,nele/2
           P(i,j) = 2.d0*C(i,k) * C(j,k) + P(i,j)
           !This will give: P(i,j) = 2* SUM[C(i,k)*C(j,k)] summing from
           !k=1 to n/2 and i,j =1 to n. This equation is from EQ 3.145 in the Szabo book
        end do
     end do
  end do
  ! print*, "Guess P"
  ! Do i = 1, N
  !    print*, (P(i,j),j=1,n)
  ! end Do
  ! open(22,file='NewDensity.txt')
  ! write(22,*) n
  ! write(22,*) P
  ! !Writing out the index size n and then the density matrix P
  ! close(21)
  ! close(22)
End Subroutine NewDensity
!
!**********************************************************************
Subroutine Sort(N,E,Vect)
  ! Sort Eigenvalues and Eigenvectors in energy order 
  implicit none
  Double precision, Parameter :: large = 1.d12
  integer :: n,i,j,savej
  integer, allocatable, dimension(:) :: Mapi, Mapj
  double precision, allocatable, dimension(:,:) :: Scr
  double precision, dimension(N,N) :: E, Vect
  double precision small
  Allocate(Scr(N,N),Mapi(N),Mapj(N))
  ! Order Eigenvalues and create map
  Mapi = 0
  Mapj = 0
  Scr = 0.0
  Do i = 1, N
     small = large
     Do j = 1, N
        If(E(j,j)< small .and. Mapj(j)==0) then
           savej = j
           small = E(j,j)
        end If
     end Do
     Mapi(i) = savej
     Mapj(savej) = 1
     ! carrying around the enrtire matrix is not very efficient and should be fixed
     Scr(i,i) = small
  end Do
  E = Scr
  ! Sort Eigenvectors
  ! print*, "Vect bef"
  ! Do i = 1, N
  !    print*, (vect(i,j),j=1,n)
  ! end Do
  Do i = 1, N
     Scr(:,i) = Vect(:,Mapi(i))
  end Do
  Vect = Scr
  ! print*, "Vect aft"
  ! Do i = 1, N
  !    print*, (vect(i,j),j=1,n)
  ! end Do
  Deallocate(Scr,Mapi,Mapj)
  !
End Subroutine Sort
!
!**********************************************************************
Subroutine diis
  !This program has two major purposes - (1) checking convergence and
  !(2) DIIS extrapolation. Specifically:
  !
  ! - after the fourth iteration, it tests for convergence based on the error between
  !   the current iteration and the previous one
  ! - if convergence has not been achieved, the program undergoes DIIS extrapolation
  !    and tests for convergence again
  ! - if either test concludes convergence, it immediately flags the
  !   main program to cease iterations
  ! - if neither test yields convergence, a new density matrix guess is written based
  !   on the DIIS coefficients, and the main program reiterates from there

  implicit none
  integer :: n, i, j, k, en, m, ierror1, ierror2, ierror3, ierror4
  real :: norm, detB
  character(len=3) :: canwefinallystopnow
  double precision, allocatable, dimension (:,:) :: P, P1, P2, P3, P4, B, bb, mm, invB, newP
  double precision, allocatable, dimension (:) :: e21, e32, e43, Er, Erb, eig, cvals, e21r, e32r, e43r, resid

  canwefinallystopnow = 'no'

  !This first part attempts to verify the prior existence of trials 1-4
!!!!!if the named file does not yet exist, it flags it with Iostat
  open(23,file='NewDensity1.txt',Status = 'old', Iostat=ierror1)
  open(24,file='NewDensity2.txt',Status = 'old', Iostat=ierror2)
  open(25,file='NewDensity3.txt',Status = 'old', Iostat=ierror3)
  open(26,file='NewDensity4.txt',Status = 'old', Iostat=ierror4)
  close(23)
  close(24)
  close(25)
  close(26)

  !Using those Iostat flags, the script will proceed to direct the most recent 
  !    iterated step ('NewDensity.txt') to its proper position in memory 
  IF (ierror1 > 0) THEN  
     !e.g. if the program was unable to open 'NewDensity1.txt' gracefully, it signals the
     !    script to write the current iteration as #1
     open(22,file='NewDensity.txt') 
     read(22,*) n
     allocate(P(n,n))
     read(22,*) P
     open(23,file='NewDensity1.txt') 
     write(23,*) n
     write(23,*) P    
     close(22)
     close(23)
  ELSE IF (ierror2 > 0) THEN   
     !if #2 doesn't exist, then write it
     open(22,file='NewDensity.txt') 
     read(22,*) n
     allocate(P(n,n))
     read(22,*) P
     open(24,file='NewDensity2.txt') 
     write(24,*) n
     write(24,*) P
     close(22)
     close(24)
  ELSE IF (ierror3 > 0) THEN   
     !if #3 doesnt exist.....etc.
     open(22,file='NewDensity.txt') 
     read(22,*) n
     allocate(P(n,n))
     read(22,*) P
     open(25,file='NewDensity3.txt') 
     write(25,*) n
     write(25,*) P
     close(22)
     close(25)
  ELSE IF (ierror4 > 0) THEN   
     open(22,file='NewDensity.txt') 
     read(22,*) n
     allocate(P(n,n))
     read(22,*) P
     open(26,file='NewDensity4.txt') 
     write(26,*) n
     write(26,*) P
     close(22)
     close(26)
  ELSE   
     !if all #1-4 already exist, then the newest iteration will be stored as #4, meaning
     !    the current #4 will replace #3, the current #3 replaces #2, and #2 replaces #1 
     open(24,file='NewDensity2.txt') 
     read(24,*) n
     allocate(P(n,n))
     read(24,*) P
     open(23,file='NewDensity1.txt') 
     write(23,*) n
     write(23,*) P   
     close(24)
     close(23)
!!!!!changes P1 -> P2
     open(25,file='NewDensity3.txt') 
     read(25,*) n
     read(25,*) P
     open(24,file='NewDensity2.txt') 
     write(24,*) n
     write(24,*) P   
     close(24)
     close(25)
!!!!!changes P2 -> P3
     open(26,file='NewDensity4.txt') 
     read(26,*) n
     read(26,*) P
     open(25,file='NewDensity3.txt') 
     write(25,*) n
     write(25,*) P   
     close(26)
     close(25)
!!!!!changes P3 -> P4
     open(22,file='NewDensity.txt') 
!!!!!and lastly use most recent P as P4
     read(22,*) n
     read(22,*) P
     open(26,file='NewDensity4.txt') 
     write(26,*) n
     write(26,*) P
     close(22)
     close(26)
  END IF
  !Thus the above block ensures that the DIIS procedure carries no more than has 4 
  !    iterations in its memory


  !The next block is the initial convergence test, which is not applied until at least
  !    four iterations have been calculated 
  IF (ierror4 == 0) THEN  
!!!!!if indeed there are at least 4, the latest two (#3 and #4) are extracted to 
     !    generate the error vector
     open(26,file='NewDensity4.txt')
     read(26,*) n
     allocate(P4(n,n))
     read(26,*) P4
     open(25,file='NewDensity3.txt')
     read(25,*) n
     allocate(P3(n,n))
     read(25,*) P3
     en = n * n 
     allocate(e43(en))
     k = 1
     IF (k <= en) THEN 
        do j=1,n
           do i=1,n
              e43(k) = P4(i,j) - P3(i,j)
              k = k + 1 
           end do
        end do
        close(25)
        close(26)
     END IF
!!!!!the error vector is then squared and its norm - the value by which convergence
     !    is measured - is determined
     allocate(Er(en))
     do k=1,en 
        Er(k) = e43(k) * e43(k)
     end do
     norm = SQRT(SUM(Er))
     print*, "norm1:", norm
!!!!!if the norm of the error vector is below a set threshold, the program signals
     !    convergence 
     IF (norm <= 0.000001) THEN
        canwefinallystopnow = 'yes'
     END IF
     print*, "can we?:", canwefinallystopnow
  END IF

  !If convergence has not been reached, the program proceeds with DIIS extrapolation by
  !    generating residual vectors between each of the four iteration
  IF (canwefinallystopnow /= 'yes') THEN 
     open(24,file='NewDensity2.txt')  
     read(24,*) n
     allocate(P2(n,n))
     read(24,*) P2
     allocate(e32(en))
     k = 1
     IF (k <= en) THEN
        do j=1,n
           do i=1,n
              e32(k) = P3(i,j) - P2(i,j)
              k = k + 1 
           end do
        end do
     END IF
     open(23,file='NewDensity1.txt')  
     read(23,*) n
     allocate(P1(n,n))
     read(23,*) P1
     allocate(e21(en))
     k = 1
     IF (k <= en) THEN
        do j=1,n
           do i=1,n
              e21(k) = P2(i,j) - P1(i,j)
              k = k + 1 
           end do
        end do
     END IF
     close(24)
     close(23)
!!!!!with the residuals e21-e43 defined, we now construct the B matrix that consists
     !    of their overlaps:
     allocate(Erb(en))
     allocate(B(4,4))
     m = 3
     allocate(bb(m,en))
!!!!!first the three residual vectors are packaged as an array for convenient 
     !    calculation of overlap
     do i=1,en
        bb(1,i) = e21(i)
        bb(2,i) = e32(i)
        bb(3,i) = e43(i)
     end do
     do j=1,4
        do i=1,4
           IF ((i <= m) .AND. (j <= m)) THEN    
              do k=1,en
                 Erb(k) = bb(i,k) * bb(j,k) 
              end do
              B(i,j) = SUM(Erb)
!!!!!the conditions below serve to populate the extra dimension of the B matrix, 
              !    which satisfies the requirement that the sum of the calculated coefficients do
              !    not exceed 1
           ELSE IF (i /= j) THEN
              B(i,j) = -1
           ELSE
              B(i,j) = 0
           END IF
        end do
     end do


     !With the B matrix constructed, the next step is to solve the eigenvalue equation.
     !This is accomplished by first determining the inverse of B, then applying it to the
     !    right side vector (0,0,0,-1)
!!!!!Calculating the inverse requires two steps - finding the determinant of B, and 
     !    obtaining the adjoint matrix of cofactors. Since the determinant can be defined
     !    by the terms of the cofactor matrix, the latter is calculated first:
     allocate(mm(4,4))
     mm(1,1) = (B(2,2) * B(3,3) * B(4,4)) - (B(2,2) * B(3,4) * B(4,3)) + &
          &(B(2,3) * B(3,4) * B(4,2)) - (B(2,3) * B(3,2) * B(4,4)) &
          &+ (B(2,4) * B(3,2) * B(4,3)) - (B(2,4) * B(3,3) * B(4,2)) 
     mm(1,2) = (B(2,1)*B(3,4)*B(4,3)) - (B(2,1)*B(3,3)*B(4,4)) + &
          &(B(2,3)*B(3,1)*B(4,4)) - (B(2,3)*B(3,4)*B(4,1)) + &
          &(B(2,4)*B(3,3)*B(4,1)) - (B(2,4)*B(3,1)*B(4,3)) 
     mm(1,3) = (B(2,1)*B(3,2)*B(4,4)) - (B(2,1)*B(3,4)*B(4,2)) + &
          &(B(2,2)*B(3,4)*B(4,1)) - (B(2,2)*B(3,1)*B(4,4)) + &
          &(B(2,4)*B(3,1)*B(4,2)) - (B(2,4)*B(3,2)*B(4,1)) 
     mm(1,4) = (B(2,1)*B(3,3)*B(4,2)) - (B(2,1)*B(3,2)*B(4,3)) + &
          &(B(2,2)*B(3,1)*B(4,3)) - (B(2,2)*B(3,3)*B(4,1)) + &
          &(B(2,3)*B(3,2)*B(4,1)) - (B(2,3)*B(3,1)*B(4,2)) 
     mm(2,1) = (B(1,2)*B(3,4)*B(4,3)) - (B(1,2)*B(3,3)*B(4,4)) + &
          &(B(1,3)*B(3,2)*B(4,4)) - (B(1,3)*B(3,4)*B(4,2)) &
          &+ (B(1,4)*B(3,3)*B(4,2)) - (B(1,4)*B(3,2)*B(4,3)) 
     mm(2,2) = (B(1,1)*B(3,3)*B(4,4)) - (B(1,1)*B(3,4)*B(4,3)) + &
          &(B(1,3)*B(3,4)*B(4,1)) - (B(1,3)*B(3,1)*B(4,4)) &
          &+ (B(1,4)*B(3,1)*B(4,3)) - (B(1,4)*B(3,3)*B(4,1)) 
     mm(2,3) = (B(1,1)*B(3,4)*B(4,2)) - (B(1,1)*B(3,2)*B(4,4)) + &
          &(B(1,2)*B(3,1)*B(4,4)) - (B(1,2)*B(3,4)*B(4,1)) &
          &+ (B(1,4)*B(3,2)*B(4,1)) - (B(1,4)*B(3,1)*B(4,2)) 
     mm(2,4) = (B(1,1)*B(3,2)*B(4,3)) - (B(1,1)*B(3,3)*B(4,2)) + &
          &(B(1,2)*B(3,3)*B(4,1)) - (B(1,2)*B(3,1)*B(4,3)) &
          &+ (B(1,3)*B(3,1)*B(4,2)) - (B(1,3)*B(3,2)*B(4,1)) 
     mm(3,1) = (B(1,2)*B(2,3)*B(4,4)) - (B(1,2)*B(2,4)*B(4,3)) + &
          &(B(1,3)*B(2,4)*B(4,2)) - (B(1,3)*B(2,2)*B(4,4)) &
          &+ (B(1,4)*B(2,2)*B(4,3)) - (B(1,4)*B(2,3)*B(4,2)) 
     mm(3,2) = (B(1,1)*B(2,4)*B(4,3)) - (B(1,1)*B(2,3)*B(4,4)) + &
          &(B(1,3)*B(2,1)*B(4,4)) - (B(1,3)*B(2,4)*B(4,1)) &
          &+ (B(1,4)*B(2,3)*B(4,1)) - (B(1,4)*B(2,1)*B(4,3)) 
     mm(3,3) = (B(1,1)*B(2,2)*B(4,4)) - (B(1,1)*B(2,4)*B(4,2)) + &
          &(B(1,2)*B(2,4)*B(4,1)) - (B(1,2)*B(2,1)*B(4,4)) &
          &+ (B(1,4)*B(2,1)*B(4,2)) - (B(1,4)*B(2,2)*B(4,1)) 
     mm(3,4) = (B(1,1)*B(2,3)*B(4,2)) - (B(1,1)*B(2,2)*B(4,3)) + &
          &(B(1,2)*B(2,1)*B(4,3)) - (B(1,2)*B(2,3)*B(4,1)) &
          &+ (B(1,3)*B(2,2)*B(4,1)) - (B(1,3)*B(2,1)*B(4,2)) 
     mm(4,1) = (B(1,2)*B(2,4)*B(3,3)) - (B(1,2)*B(2,3)*B(3,4)) + &
          &(B(1,3)*B(2,2)*B(3,4)) - (B(1,3)*B(2,4)*B(3,2)) &
          &+ (B(1,4)*B(2,3)*B(3,2)) - (B(1,4)*B(2,2)*B(3,3)) 
     mm(4,2) = (B(1,1)*B(2,3)*B(3,4)) - (B(1,1)*B(2,4)*B(3,3)) + &
          &(B(1,3)*B(2,4)*B(3,1)) - (B(1,3)*B(2,1)*B(3,4)) &
          &+ (B(1,4)*B(2,1)*B(3,3)) - (B(1,4)*B(2,3)*B(3,1)) 
     mm(4,3) = (B(1,1)*B(2,4)*B(3,2)) - (B(1,1)*B(2,2)*B(3,4)) + &
          &(B(1,2)*B(2,1)*B(3,4)) - (B(1,2)*B(2,4)*B(3,1)) &
          &+ (B(1,4)*B(2,2)*B(3,1)) - (B(1,4)*B(2,1)*B(3,2)) 
     mm(4,4) = (B(1,1)*B(2,2)*B(3,3)) - (B(1,1)*B(2,3)*B(3,2)) + &
          &(B(1,2)*B(2,3)*B(3,1)) - (B(1,2)*B(2,1)*B(3,3)) &
          &+ (B(1,3)*B(2,1)*B(3,2)) - (B(1,3)*B(2,2)*B(3,1)) 
     !Now the determinant can be far more simply defined:
     detB = (B(1,1)*mm(1,1)) + (B(2,1)*mm(1,2)) + (B(3,1)*mm(1,3)) + (B(4,1)*mm(1,4))
     allocate(invB(4,4))
     invB = (1/(detB))*mm
     !With the inverse matrix of B calculated, we can now solve for the linear 
     !    combination vector. First we create the eigenvalue (0,0,0,-1), and then a 
     !    matrix multiplication is carried out:
     allocate(eig(4))    
     do k=1,m
        eig(k) = (0)
     end do
     eig(4) = (-1)
     allocate(cvals(4))
     do k=1,4
        cvals(k) = (invB(4,k) * eig(4)) 
     end do

     !The next step is to construct/redefine a new residuum vector based on a linear 
     !    combination of weighted residual vectors
     allocate(resid(en))
     allocate(e43r(en))
     allocate(e32r(en))
     allocate(e21r(en))
     do k=1,en
        e43r(k) = (cvals(3))*(e43(k))
        e32r(k) = (cvals(2))*(e32(k))
        e21r(k) = (cvals(1))*(e21(k))
     end do
     do k=1,en
        resid(k) = e21r(k) + e32r(k) + e43r(k)  
     end do

     !Now we check the new residuum vector for convergence and, if it does, the script
     !    signals the main program
     do k=1,en 
        Er(k) = resid(k) * resid(k)
     end do
     norm = SQRT(SUM(Er))
     print*, "norm2:", norm
     IF (norm <= 0.000001) THEN
        canwefinallystopnow = 'yes'
        !If convergence has still not been achieved, a new density
        !    matrix is constructed as a linear combination of the
        !    previously iterated matrices weighted by the DIIS
        !    coefficients, and subsequently exports them.
     ELSE
        allocate(newP(n,n))
        do j=1,n
           do k=1,n
              newP(k,j) = (P2(k,j)*cvals(1)) + (P3(k,j)*cvals(2)) + (P4(k,j)*cvals(3))   
           end do
        end do
        open(28,file='NewDensity.txt')
        write(28,*) n
        write(28,*) newP
        close(28) 
     END IF
     print*, "can we now?:", canwefinallystopnow
  END IF
  !
End Subroutine diis
!
!**********************************************************************
Function Trace(n,A,B)
  ! This routine evaluate the trace of two square matrices of dimension n
  implicit none
  integer i,j,n
  double precision, dimension (n,n) :: A, B
  double precision t, Trace
  t = 0.d0  
  Do i = 1, n
     Do j = 1, n
        t = t + A(i,j)*B(j,i)
     end Do
  end Do
  Trace = t
End Function trace
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************

Program KUHF
  Use KUHF_Integrals
  !
  ! This program performs a Hartree-Fock (HF) calculation.
  !
  ! NB : number of basis functions
  ! Vnn : nuc-nuc repulsion energy
  implicit none
  integer NB, Nele, iSCF, maxcycle 
  integer i,j
  Integer, Parameter:: NBMax = 100
  double precision Vnn, DeltaE, Energy, OldEn, DeltaP, ThrE, ThrP, Trace
  double precision, allocatable, dimension(:,:) :: Core, S, Fock, P, X, G, C, OldP
1  Format(F10.6,1X)
10 Format(1x,'SCF Cycle ',I3,3x,'Total Energy=',F20.10,3x,'Delta-E=',F20.10)
11 Format(1x,'Final Energy=',F20.10,3x,'Delta-E=',F20.10)
20 Format(1x,'Nuclear Repulsion Energy=',F20.10)
  open(10, file='kuhf.txt')
  !
  ThrE = 1.D-9
  ThrP = 1.D-7
  ! Read-in geometry and compute integrals in AO basis
  Call Integrals(NB,Nele,Vnn,S,Core)
  If(NB > NBMax) then
     print*, "Basis set size exceeds Max"
     stop
  end If  
  allocate (Fock(NB,NB))
  allocate (P(NB,NB))
  allocate (X(NB,NB))
  allocate (G(NB,NB))
  allocate (C(NB,NB))
  allocate (OldP(NB,NB))
  !  print core
  !  print*, "Core"
  ! Do i = 1, NB
  !    print*, (Core(i,j),j=1,nb)
  ! end Do
  Write(*,20) Vnn
  ! Diagnalize the overal matrix and compute the tranformation matrix
  ! for the canonical orthogonalization
  ! print*, "S1"
  ! Do i = 1, NB
  !    print*, (S(i,j),j=1,nb)
  ! end Do
  Call SDiag(NB,S,X)
  ! print*, "S2"
  ! Do i = 1, NB
  !    print*, (S(i,j),j=1,nb)
  ! end Do
  Call TransformationMatrix(NB,S,X)
  ! print*, "S3"
  ! Do i = 1, NB
  !    print*, (S(i,j),j=1,nb)
  ! end Do
  ! print*, "X"
  ! Do i = 1, NB
  !    print*, (X(i,j),j=1,nb)
  ! end Do
  Deallocate(S)
  ! Produce first density guess fro core Hamiltonian
  Fock = Core
  Call transformedF(NB,Fock,X)
  ! print*, "F-Core0"
  ! Do i = 1, NB
  !    print*, (Fock(i,j),j=1,nb)
  ! end Do
  Call SDiag(NB,Fock,C)
!  Call FDiag(NB,Fock,C)
  ! print*, "F-Core"
  ! Do i = 1, NB
  !    print*, (Fock(i,j),j=1,nb)
  ! end Do
  ! print*, "Guess Cp"
  ! Do i = 1, NB
  !    print*, (C(i,j),j=1,nb)
  ! end Do
  Call Cmatrix(NB,X,C)
!   print*, "Guess C"
!   Do i = 1, NB
!      print*, (C(i,j),j=1,nele/2)
! !     print*, (C(i,j),j=1,nb)
!   end Do
  Call NewDensity(NB,Nele,C,P)
  Call Gmatrix(NB,G,P)
  Fock = G + Core
!   print*, "Guess Density"
!   Do i = 1, NB
!      print*, P(i,1)
! !     print*, (P(i,j),j=1,nb)
!   end Do
  Energy = Trace(NB,P,Core) 
  print*,'E(core)',energy,0.5*energy
  !
  ! Now do the actual loop
  iSCF = 1
  maxcycle = 50
  OldEn = 0.d0
  DeltaE = 1.d0
  DeltaP = 1.d0
  OldP = P
  DeltaP = Trace(NB,P,P)
!  write(*,*)'p guess',deltap
  Do While (Abs(DeltaE) >= ThrE .and. DeltaP >= ThrP .and. iSCF <= maxcycle ) 
     !     print*, "enters",iscf
     ! print core
     ! print*, "Density"
     ! Do i = 1, NB
     !    print*, (P(i,j),j=1,nb)
     ! end Do
     ! print*, "G"
     ! Do i = 1, NB
     !    print*, (G(i,j),j=1,nb)
     ! end Do     
     ! print*, "F"
     ! Do i = 1, NB
     !    print*, (Fock(i,j),j=1,nb)
     ! end Do     
     Energy = 0.5d0 *Trace(NB,P,G) + Trace(NB,P,Core) + Vnn
     DeltaE = Energy - OldEn
     Write(*,10) ,iSCF,Energy,DeltaE
!     stop
     OldEn = Energy
     Call transformedF(NB,Fock,X)
     Call SDiag(NB,Fock,C)
!     Call FDiag(NB,Fock,C)
     ! Sort Eigenvalues and eigenvectors in increasing order
     Call Sort(NB,Fock,C)
     Call Cmatrix(NB,X,C)
     Call NewDensity(NB,Nele,C,P)
     Call Gmatrix(NB,G,P)
     Fock = G + Core
     OldP = P - OldP
     DeltaP = Trace(NB,OldP,OldP)
!     write(*,*)'deltap',deltap
     If(DeltaP<0.d0) then
        print*,"Something is wrong with P. DeltaP=",DeltaP
        stop
     end If
     DeltaP = sqrt(deltap)
     OldP = P
     iSCF = iSCF + 1
  end Do
  ! Print final energy
  Energy = 0.5d0 *Trace(NB,P,G) + Trace(NB,P,Core) + Vnn
!  print*,Trace(NB,P,G),Trace(NB,P,Core),Vnn
  DeltaE = Energy - OldEn
  Write(*,11) ,Energy,DeltaE
  Deallocate(Fock,P,X,G,C,OldP)
  Deallocate(Core)
  close(10)
  !
End Program KUHF
