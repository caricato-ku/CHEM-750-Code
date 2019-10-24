!  This is bootstrap unrestricted CCSD program written by the Spring
!  2015 class of Advanced Quantum Mechanics (CHEM 850) at the
!  University of Kansas: Tal Aharon, Matthew Barclay, Sijin Ren, and
!  Marco Caricato.
!
!  This program is designed to accept previously calculated
!  Hartree-Fock two-electron integrals and apply the coupled-cluster
!  post-SCF method.
!
!  This program will be carried out in a number of steps:
!    (1) Reading in the previously calculated integrals
!    (2) Calculating an initial guess of the excitation Amplitudes
!          defined solely as the corresponding two-electron integrals
!    (3) Iteratively solving for the excitation Amplitudes and
!    incorporating them into the final coupled-cluster energy
!
!  This first block creates a Module for Global variables, allowing us
!  to separate our code into a number of Subroutines that each
!  calculate various Terms of the coupled-cluster Amplitude equation.
!
! Input file:
! Method = choice of method (CCD or CCSD)
! EpsI = occupied orbital energies
! EpsA = virtual orbital energies
! abcd = <ab||cd> integrals (in that index order)
! iabc = <ia||bc> integrals (in that index order)
! iajb = <ia||jb> integrals (in that index order)
! ijab = <ij||ab> integrals (in that index order)
! ijka = <ij||ka> integrals (in that index order)
! ijkl = <ij||kl> integrals (in that index order)
!
! The name of the files is generic. If no file name is specified, the
! default names are: EpsA, EpsI, Iabcd, Iiabc, Iiajb, Iijab, Iijka, Iijkl.
!
Module Global
   Double Precision, Allocatable, Dimension(:,:,:,:) :: abcd, iabc, iajb, ijab, ijka, ijkl, DubAmpl, Dijab, AmplDiff
   Double Precision, Allocatable, Dimension (:,:) :: SingAmpl, Dia
   Double Precision, Allocatable, Dimension (:) :: EpsI, EpsA
End Module Global
!
Program Overall
   Implicit None
   Character(300) :: Method, Junk
   Integer IError
   Open(18, File='Input', Status='old', IOStat = IError)
   If (IError .ne. 0) Then
     Print*, "No Input file found"
     Call Exit() 
   End If
   Read(18,*) Junk, Junk, Method
   Close(18)
   If (Method == "CCSD") Then
     Call CCSD
   else If (Method == "CCD") Then
     Call CCD
   else
      Write(6,*) "Method Not Recognized"
   End If
End Program Overall
!
!!This next block houses the actual program, which Calls the various Subroutines and calculates the Amplitudes 
!    and Energy.
Subroutine CCSD
   Use Global
   Implicit Integer (a-n)
   Character(300) :: FileName, Junk
   Double Precision, Allocatable, Dimension(:,:,:,:) :: UpdDubAmpl
   Double Precision, Allocatable, Dimension(:,:) :: UpdSingAmpl
   Double Precision ConvA, ConvE, Energy, PrevEnergy, Term1, Term2, Term3, Term4, Term5, Term6,&
&Term7, Term8, Term9, Term10, Retval, Ret, Temp, Temp2, Term11, Term12, Term13, Term14,&
&TermA, TermB, TermC, TermD, TermE, TermF, Temp3
   Integer NO, NV, aaa, Counter
   Common NO, NV
!
!!First, the program reads the previously calculated values based on occupied (i,j,k,l...m,n) orbitals and virtual (a,b,c,d...e,f) orbitals:
   Open(40, File='Input')
   Read(40,*)
   Read(40,*,IOStat=IError3) Junk, Junk, FileName
!EpsI  !default if it doesn't find the file is EpsI  
   Open(18, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(18,*) NO
     Allocate(EpsI(NO))
     Do i = 1, NO
        Read(18,*) EpsI(i)
     End Do
     IError2 = 0
   End If 
   Close(18)
   If (IError .ne. 0) Then
     Open(18, File='EpsI', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(18,*) NO
       Allocate(EpsI(NO))
       Do i = 1, NO
          Read(18,*) EpsI(i)
       End Do
     End If
     Close(18)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "EpsI not found"
     Return
   End If
!!!!!where we have Read the orbital energies (EpsI) of the Occupied orbitals
   Read(40,*,IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!EpsA  !default if it doesn't find the file is EpsA  
   Open(9, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(9,*) NV
     Allocate(EpsA(NV))
     Do i = 1, NV
        Read(9,*) EpsA(i)
     End Do
     IError2=0
   End If 
   Close(9)
   If (IError .ne. 0) Then
     Open(9, File='EpsA', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(9,*) NV
       Allocate(EpsA(NV))
       Do i = 1, NV
          Read(9,*) EpsA(i)
       End Do
     End If
     Close(9)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "EpsA not found"
     Return
   End If
!!!!!similarly, we have Read the energies of the Virtual orbitals
!!!!!Next, the two-electron integrals (of form <wx||yz>) are Read:
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iabcd  !default if it doesn't find the file is Iabcd 
   Open(11, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(11,*) NO, NV
     Allocate(abcd(NV,NV,NV,NV))
     Do d = 1, NV
        Do c = 1, NV
           Do b = 1, NV
              Do a = 1, NV
                 Read(11,*) abcd(a,b,c,d)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(11)
   If (IError .ne. 0) Then
     Open(11, File='Iabcd', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(11,*) NO, NV
       Allocate(abcd(NV,NV,NV,NV))
       Do d = 1, NV
          Do c = 1, NV
             Do b = 1, NV
                Do a = 1, NV
                   Read(11,*) abcd(a,b,c,d)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(11)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ab||cd> not found"
     Return
   End If
!!!!!here we have constructed the Global 4-Dimensional matrix of integrals whose indices run over 
!         only virtual orbitals (e.g. a,b,c,d)
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iiabc !default if it doesn't find the file is Iiabc 
   Open(12, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(12,*) NO, NV
     Allocate(iabc(NO,NV,NV,NV))
     Do c = 1, NV
        Do b = 1, NV
           Do a = 1, NV
              Do i = 1, NO
                 Read(12,*) iabc(i,a,b,c)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(12)
   If (IError .ne. 0) Then
     Open(12, File='Iiabc', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(12,*) NO, NV
       Allocate(iabc(NO,NV,NV,NV))
       Do c = 1, NV
          Do b = 1, NV
             Do a = 1, NV
                Do i = 1, NO
                   Read(12,*) iabc(i,a,b,c)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(12)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ia||bc> not found"
     Return
   End If
!!!!!here we have constructed a similar matrix of integrals involving three virtual orbitals and one occupied orbital 
!
!!!!!similarly, the Global variable of two occupied, two virtual integrals in the form <ia||jb>:
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iiajb !default if it doesn't find the file is Iiajb  
   Open(13, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(13,*) NO, NV
     Allocate(iajb(NO,NV,NO,NV))
     Do b = 1, NV
        Do j = 1, NO
           Do a = 1, NV
              Do i = 1, NO
                 Read(13,*) iajb(i,a,j,b)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(13)
   If (IError .ne. 0) Then
     Open(13, File='Iiajb', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(13,*) NO, NV
       Allocate(iajb(NO,NV,NO,NV))
       Do b = 1, NV
          Do j = 1, NO
             Do a = 1, NV
                Do i = 1, NO
                   Read(13,*) iajb(i,a,j,b)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(13)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ia||jb> not found"
     Return
   End If
!!!!!and the same for the form <ij||ab>:
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iijab !default if it doesn't find the file is Iijab  
   Open(14, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(14,*) NO, NV
     Allocate(ijab(NO,NO,NV,NV))
     Do b = 1, NV
        Do a = 1, NV
           Do j = 1, NO
              Do i = 1, NO
                 Read(14,*) ijab(i,j,a,b)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(14)
   If (IError .ne. 0) Then
     Open(14, File='Iijab', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(14,*) NO, NV
       Allocate(ijab(NO,NO,NV,NV))
       Do b = 1, NV
          Do a = 1, NV
             Do j = 1, NO
                Do i = 1, NO
                   Read(14,*) ijab(i,j,a,b)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(14)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ij||ab> not found"
     Return
   End If
!!!!!and for three occupied, one virtual (<ij||ka>):
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iijka !default if it doesn't find the file is Iijka  
   Open(15, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(15,*) NO, NV
     Allocate(ijka(NO,NO,NO,NV))
     Do a = 1, NV
        Do k = 1, NO
           Do j = 1, NO
              Do i = 1, NO
                 Read(15,*) ijka(i,j,k,a)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(15)
   If (IError .ne. 0) Then
     Open(15, File='Iijka', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(15,*) NO, NV
       Allocate(ijka(NO,NO,NO,NV))
       Do a = 1, NV
          Do k = 1, NO
             Do j = 1, NO
                Do i = 1, NO
                   Read(15,*) ijka(i,j,k,a)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(15)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ij||ka> not found"
     Return
   End If
!!!!!and for integrals of only occupied orbitals:
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iijkl !default if it doesn't find the file is Iijkl  
   Open(16, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(16,*) NO, NV
     Allocate(ijkl(NO,NO,NO,NO))
     Do l = 1, NO
        Do k = 1, NO
           Do j = 1, NO
              Do i = 1, NO
                 Read(16,*) ijkl(i,j,k,l)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(16)
   If (IError .ne. 0) Then
     Open(16, File='Iijkl', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(16,*) NO, NV
       Allocate(ijkl(NO,NO,NO,NO))
       Do l = 1, NO
          Do k = 1, NO
             Do j = 1, NO
                Do i = 1, NO
                   Read(16,*) ijkl(i,j,k,l)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(16)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ij||kl> not found"
     Return
   End If
   Close(40)
!
!!The next block of the main program constructs the Amplitudes and
!!ultimately deTermines the Energy.
!!!!!First we Allocate memory to the 4-Dimensional matrices for
!!!!!Amplitudes and the denominator Terms 'Dijab' and 'Dia'
!         which (we see below) is the sum of occupied orbital energies
!         minus the virtual orbitals they will be excited to in the
!         corresponding cluster excitation
   Allocate(UpdDubAmpl(NO,NO,NV,NV))
   Allocate(DubAmpl(NO,NO,NV,NV))
   Allocate(Dijab(NO,NO,NV,NV))
   Allocate(AmplDiff(NO,NO,NV,NV))
   Allocate(Dia(NO,NV))
   Allocate(SingAmpl(NO,NV))
   Allocate(UpdSingAmpl(NO,NV))
!   
!!!!!Thus we calculate the Dijab and Dia matrices:
   Do b = 1, NV
      Do a = 1, NV
         Do j = 1, NO
            Do i = 1, NO
               Dijab(i,j,a,b) = EpsI(i) + EpsI(j) - EpsA(a) - EpsA(b)
               Dia(i,a) = EpsI(i) - EpsA(a)
            End Do
         End Do
      End Do
   End Do
!
!!Here we perform our initial guess of the excitation Amplitudes based
!!solely on the two-electron integrals from
!    the SCF calculation: 
   Do b = 1, NV
      Do a = 1, NV
         Do j = 1, NO
            SingAmpl(j,a) = 0
            Do i = 1, NO
               DubAmpl(i,j,a,b) = ijab(i,j,a,b)/Dijab(i,j,a,b)
            End Do
         End Do
      End Do
   End Do
!!!!!We Then perform an initial Energy calculation based on these
!!!!!Amplitude guesses, giving the MP2 energy:
   Do b = 1, NV
      Do a = 1, NV
         Do j = 1, NO
            Do i = 1, NO
               Call Tau(DubAmpl(i,j,a,b),i,j,a,b, Retval)
               Energy = Energy + .25*Retval*ijab(i,j,a,b)
!!!!!where 'Tau' is a Subroutine function of the excitation Amplitudes
            End Do
         End Do
      End Do
   End Do
   print*, "MP2 Energy", Energy
!
!!The program now proceeds to perform the actual iterations.
!!!!!Each of the numbered "Terms" (Term1, Term2, Term3...) is one of
!!!!!the parts of the Double-excitation cluster operator (T2)
!!!!!equation, as it is broken Down and explicitly defined before
!!!!!being recombined at the End.
   Energy = 0
   PrevEnergy = 0
   aaa = 0
   Counter = 0
!!!!!thus we have provided an initial guess of zero for each of our
!!!!!following variables, allowing for ease of summation
!
!
!!!!!The following block of code is the main structure of the
!!!!!iterative loop,including the test of
!         and criteria for Convergence.
   Do
      If (Counter == 100) Then
         print*, "Convergence not achieved. Current energy after 100 cycles is ", Energy
         Exit
      End If
      If (aaa == 1) Then !Convergence
         Exit
      End If
!!!!!where we have initially determined that If the program Does not Converge within 100 iterations, it will cease.
      Do b = 1, NV
         Do a = 1, NV
            Do j = 1, NO
               Do i = 1, NO
                  Term1 = 0
                  Term2 = 0
                  Term3 = 0
                  Term4 = 0
                  Term5 = 0
                  Term6 = 0
                  Term7 = 0
                  Term8 = 0
                  Term9 = 0
                  Term10 = 0
                  Term11 = 0
                  Term12 = 0
                  Term13 = 0
                  Term14 = 0
!!!!!here we have constructed the iterative loops over all indices of occupied (i,j) and virtual (a,b) orbitals 
!         prior to actually defining the individual Terms
!!Terms 1 and 2 consist of two permutations of the virtual orbitals (ab) within a sum over a single virtual orbital
!    (e) of Terms of the "F" matrix [defined as a Subroutine] whose corresponding indices are virtual orbitals: 
                  Do e = 1, NV
                     Temp = 0
                     Temp2 = 0
                     Do m = 1, NO
                        Call FBoth(m, e, Retval)
                        Temp = Temp - ((.5)*SingAmpl(m,b)*Retval)
                        Temp2 = Temp2 - ((.5) * SingAmpl(m,a)*Retval)
!mc                        Temp2 = Temp2 + ((.5) * SingAmpl(m,a)*Retval)
                     End Do
                     Call FVirt(b,e,Retval)
                     Term1 = Term1 + DubAmpl(i,j,a,e) * Retval + DubAmpl(i,j,a,e)*Temp
                     Call FVirt(a,e,Retval)
                     Term2 = Term2 + DubAmpl(i,j,b,e) * Retval + DubAmpl(i,j,b,e)*Temp2
                  End Do
!!!!!thus we have summed every (*,e) element of the virtual-orbital F matrix ['FVirt'], constructing the Terms 1 and 2
!         [F(b,e) and F(a,e) respectively] alongside the corresponding excitation Amplitudes
!!Similarly, Terms 3 and 4 consist of the contributions where a single occupied orbital (m) dIffers in Terms containing
!    the occupied-orbital F matrix ['FOcc']; each Term is a permutation of the occupied orbitals (i,j): 
                  Do m = 1, NO
                     Temp = 0
                     Temp2 = 0
                     Do e = 1, NV
                        Call FBoth(m, e, Retval)
                        Temp = Temp + ((.5)*SingAmpl(j,e)*Retval)
                        Temp2 = Temp2 + ((.5) * SingAmpl(i,e)*Retval)
!mc                        Temp2 = Temp2 - ((.5) * SingAmpl(i,e)*Retval)
                     End Do
                     Call FOcc(m, j, Retval)
                     Term3 = Term3 + DubAmpl(i,m,a,b) * Retval + DubAmpl(i,m,a,b) * Temp
                     Call FOcc(m, i, Retval)
                     Term4 = Term4 + DubAmpl(j,m,a,b) * Retval + DubAmpl(j,m,a,b) * Temp2
                  End Do
!!Terms 5 and 6 are simple sums over two occupied (m,n) and two virtual (e,f) orbitals respectively, and Call upon the
!    corresponding W matrices and excitation Amplitudes: 
                  Do n = 1, NO
                     Do m = 1, NO
                        Call Tau(DubAmpl(m,n,a,b),m,n,a,b, Ret)
                        Call WOcc(m,n,i,j, Retval)
                        Term5 = Term5 + .5*Ret*Retval
                     End Do
                  End Do
                  Do f = 1, NV
                     Do e = 1, NV
                        Call Tau(DubAmpl(i,j,e,f),i,j,e,f, Ret)
                        Call WVirt(a,b,e,f, Retval)
                        Term6 = Term6 + .5*Ret*Retval
                     End Do
                  End Do
!!Terms 6-10 are the contributions of the mixed (two occupied, two virtual) W matrix ['WBoth'] summed over one dIffering
!    occupied and one dIffering virtual orbital; each Term is part of simultaneous permutations of (a,b) and (i,j):
                  Do e = 1, NV
                     Do m = 1, NO
                        Call WBoth(m,b,e,j, Retval)
                        Term7 = Term7 + DubAmpl(i,m,a,e)*Retval - SingAmpl(i,e)*SingAmpl(m,a)&
                        &*(-iajb(m,b,j,e))
                        Call WBoth(m,a,e,j, Retval)
                        Term8 = Term8 + DubAmpl(i,m,b,e)*Retval  - SingAmpl(i,e)*SingAmpl(m,b)&
                        &*(-iajb(m,a,j,e))
                        Call WBoth(m,b,e,i, Retval)
                        Term9 = Term9 + DubAmpl(j,m,a,e)*Retval - SingAmpl(j,e)*SingAmpl(m,a)&
                        &*(-iajb(m,b,i,e))
                        Call WBoth(m,a,e,i, Retval)
                        Term10 = Term10 + DubAmpl(j,m,b,e)*Retval - SingAmpl(j,e)*SingAmpl(m,b)&
                        &*(-iajb(m,a,i,e))
                     End Do
                  End Do
!
                  Do e = 1, NV
                     Term11 = Term11 + SingAmpl(i,e) * iabc(j,e,b,a)
                     Term12 = Term12 + SingAmpl(j,e) * iabc(i,e,b,a)
                  End Do
!
                  Do m = 1, NO
                     Term13 = Term13 + SingAmpl(m,a) * ijka(i,j,m,b)
                     Term14 = Term14 + SingAmpl(m,b) * ijka(i,j,m,a)
                  End Do
!
!!Having defined all the Terms in the Amplitude equation, we combine them all together in the 'updated' excitation
!    Amplitude:
                  UpdDubAmpl(i,j,a,b) = (ijab(i,j,a,b)+Term1-Term2-Term3&
&+Term4+Term5+Term6+Term7-Term8-Term9+Term10+Term11-Term12-Term13+Term14)/Dijab(i,j,a,b)
              End Do
            End Do
         End Do
      End Do
!Here the Terms for the T1 Amplitude are calculated, similar to those of the T2 Amplitudes
      Do a = 1, NV
         Do i = 1, NO
            TermA = 0
            TermB = 0
            TermC = 0
            TermD = 0
            TermE = 0
            TermF = 0
            Do e = 1, NV
               Call FVirt(a,e,Retval)
               TermA = TermA + SingAmpl(i,e)*Retval
            End Do
            Do e = 1, NV
               Do m = 1, NO
                  Call FBoth(m,e,Retval)
                  TermC = TermC + DubAmpl(i,m,a,e)*Retval 
               End Do
            End Do
            Do m = 1, NO
               Call FOcc(m,i,Retval)
               TermB = TermB + SingAmpl(m,a) * Retval
            End Do
            Do f = 1, NV
               Do n = 1, NO
                  TermD = TermD + SingAmpl(n,f) * iajb(n,a,i,f)
               End Do
            End Do
            Do f = 1, NV
               Do e = 1, NV
                  Do m = 1, NO
                     TermE = TermE + .5*DubAmpl(i,m,e,f)*iabc(m,a,e,f)
                  End Do
               End Do
            End Do
            Do n = 1, NO
               Do e = 1, NV
                  Do m = 1, NO
                     TermF = TermF + .5*DubAmpl(m,n,a,e)*ijka(m,n,i,e)
                  End Do
               End Do
            End Do
            UpdSingAmpl(i,a) = (TermA - TermB + TermC - TermD - TermE - TermF)/Dia(i,a)
        End Do
      End Do
      AmplDiff = UpdDubAmpl-DubAmpl      
      !print*, UpdDubAmpl(5,6,2,3), " ", DubAmpl(5,6,2,3), " ", AmplDiff(5,6,2,3)
      DubAmpl = UpdDubAmpl
      SingAmpl = UpdSingAmpl
!!!!!and thus the Previous excitation Amplitude is redefined as the updated value
!
!!Finally, the program recalculates the coupled-cluster Energy from the new Amplitudes:
      Energy = 0
      Do b = 1, NV
         Do a = 1, NV
            Do j = 1, NO
               Do i = 1, NO
                  Call Tau(DubAmpl(i,j,a,b),i,j,a,b, Retval)
                  Energy = Energy + .25*Retval*ijab(i,j,a,b)
               End Do
            End Do
         End Do
      End Do
!
!!and checks for Convergence based on both the Energy dIfference (between current Energy and Previous Energy), as well
!    as the dIfference in the root-mean squared (RMS) value of the excitation Amplitudes:
      Call AmplConv(ConvA)
      ConvE = Energy - PrevEnergy
      ConvE = sqrt(ConvE**2)
      If ((ConvA <= 1.0D-5) .and. (ConvE <= 1.0D-7)) Then
         aaa = 1
      End If
!!!!!If Convergence has been achieved, the program signals the End of the iterative loop and prints the final Energy;
!         otherwise, it cycles back through the next iteration
      PrevEnergy = Energy
      Counter = Counter + 1
      Print*, "The energy of iteration ", Counter, "is", Energy
   End Do
   print*, "The final Energy is: ", Energy
End Subroutine CCSD
!
!
!
!!The following blocks are each of the Subroutines/Terms incorporated into the overarching Amplitude equation.
!
!
!!The Tau and "TauTilde" Subroutines calculate the values of the variables of their namesake, which are basiCally sums
!    of the operators leading to Double-excitation. 
!!!!!As this program's current function is a CCD calculation, the current values of Tau/Tautilde are solely that of the
!         corresponding Double-excitation Amplitude.
Subroutine Tau(Tijab, i, j, a, b, Retval)
   Use Global
   Implicit None
   Double Precision, Intent(in) :: Tijab
   Integer, Intent(in) :: i, a, j, b
   Double Precision, Intent(out) :: Retval
   Double Precision SingAmpDIff
   SingAmpDIff = SingAmpl(i,a)*SingAmpl(j,b) - SingAmpl(i,b)*SingAmpl(j,a)
   Retval = Tijab + SingAmpDIff
End Subroutine Tau
!
Subroutine TauTilde(Tijab, i, j, a, b, Retval)
   Use Global
   Implicit None
   Double Precision, Intent(in) :: Tijab
   Integer, Intent(in) :: i, a, j, b
   Double Precision, Intent(out) :: Retval
   Double Precision SingAmpDIff
   SingAmpDIff = (SingAmpl(i,a)*SingAmpl(j,b) - SingAmpl(i,b)*SingAmpl(j,a))/2.0
   Retval = Tijab + SingAmpDIff
End Subroutine TauTilde
!
!!The F matrices similarly correspond to their namesake and are defined by the nature of their indices - e.g. 'FVirt' 
!         is the F matrix defined only through the virtual orbitals of the excitation.
Subroutine FVirt(b, e, Retval)
   Use Global
   Implicit Integer (f-o)
   Integer, Intent(in) :: b, e
   Double Precision, Intent(out) :: Retval
   Double Precision SumMNF, TaoT, SumMF
   Integer NO, NV
   Common NO, NV
   SumMNF = 0
   SumMF = 0
   Do f = 1, NV
      Do m = 1, NO
         SumMF = SumMF + SingAmpl(m,f) * iabc(m,b,f,e)
      End Do
   End Do   
   Do f = 1, NV
      Do n = 1, NO
         Do m = 1, NO
            Call TauTilde(DubAmpl(m,n,b,f),m,n,b,f, TaoT)
            SumMNF = SumMNF + TaoT*ijab(m,n,e,f)
         End Do
      End Do
   End Do
   Retval = SumMNF * (-.5) + SumMF 
End Subroutine FVirt
!
Subroutine FOcc(m, j, Retval)
   Use Global
   Implicit None
   Integer, Intent(in) :: j, m
   Double Precision, Intent(out) :: Retval
   Double Precision SumNEF, TaoT, SumNE
   Integer n, e, f, NO, NV
   Common NO, NV
   SumNE = 0
   SumNEF = 0
   Do e = 1, NV
      Do n = 1, NO
         SumNE = SumNE + SingAmpl(n,e) * ijka(m,n,j,e)
      End Do
   End Do   
   Do f = 1, NV
      Do e = 1, NV
         Do n = 1, NO
            Call TauTilde(DubAmpl(j,n,e,f),j,n,e,f, TaoT)
            SumNEF = SumNEF + TaoT*ijab(m,n,e,f)
         End Do
      End Do
   End Do
   Retval = SumNEF * .5 + SumNE
End Subroutine FOcc
!
Subroutine FBoth(m, e, Retval)
   Use Global
   Implicit None
   Integer, Intent(in) :: e, m
   Double Precision, Intent(out) :: Retval
   Double Precision SumNF
   Integer n, f, NO, NV
   Common NO, NV
   SumNF = 0
   Do f = 1, NV
      Do n = 1, NO
         SumNF = SumNF + SingAmpl(n,f) * ijab(m,n,e,f)
      End Do
   End Do
   Retval = SumNF
End Subroutine FBoth
!
!
!!The W matrices similarly correspond to their namesake and are defined by the nature of their indices - e.g. 'WVirt' 
!         is the W matrix defined only through the virtual orbitals of the excitation.
Subroutine WOcc(m,n,i,j,Retval)
   Use Global
   Implicit None
   Integer, Intent(in) :: m, n, i, j
   Double Precision, Intent(out) :: Retval
   Double Precision SumEF, Tao, SumE, SumEPerm
   Integer e, f, NO, NV
   Common NO, NV
   SumEF = 0
   SumE = 0
   SumEPerm = 0
   Do e = 1, NV
      SumE = SumE + SingAmpl(j,e) * (ijka(m,n,i,e)) 
      SumEPerm = SumEPerm - (SingAmpl(i,e) * (ijka(m,n,j,e)))
   End Do
   Do f = 1, NV
      Do e = 1, NV
         Call Tau(DubAmpl(i,j,e,f),i,j,e,f, Tao)
         SumEF = SumEF +.25*Tao*ijab(m,n,e,f)
      End Do
   End Do
   Retval = SumEF + ijkl(m,n,i,j) + SumE + SumEPerm
End Subroutine WOcc
!
Subroutine WVirt(a,b,e,f,Retval)
   Use Global
   Implicit None
   Integer, Intent(in) :: a, b, e, f
   Double Precision, Intent(out) :: Retval
   Double Precision SumMN, Tao, SumM, SumMPerm
   Integer m, n, NO, NV
   Common NO, NV
   SumMN = 0
   SumM = 0
   SumMPerm = 0
   Do m = 1, NO
      SumM = SumM - SingAmpl(m,b) * (-iabc(m,a,e,f))
      SumMPerm = SumMPerm + (SingAmpl(m,a) * (-iabc(m,b,e,f)))
   End Do
   Do n = 1, NO
      Do m = 1, NO
         Call Tau(DubAmpl(m,n,a,b),m,n,a,b, Tao)
         SumMN = SumMN +.25*Tao*ijab(m,n,e,f)
      End Do
   End Do
   Retval = SumMN + abcd(a,b,e,f) + SumM + SumMPerm
End Subroutine WVirt
!
Subroutine WBoth(m,b,e,j,Retval)
   Use Global
   Implicit None
   Integer, Intent(in) :: m, b, e, j
   Double Precision, Intent(out) :: Retval
   Double Precision SumNF, AmpTerm, SumN, SumF
   Integer n, f, NO, NV
   Common NO, NV
   SumNF = 0
   SumN = 0
   SumF = 0
   Do f = 1, NV
      SumF = SumF + SingAmpl(j,f) * iabc(m,b,e,f)
   End Do
   Do n = 1, NO
      SumN = SumN - (SingAmpl(n,b) * ijka(n,m,j,e))
   End Do
   Do f = 1, NV
      Do n = 1, NO
         AmpTerm = (.5)*DubAmpl(j,n,f,b) + SingAmpl(j,f)*SingAmpl(n,b)
         SumNF = SumNF - AmpTerm*ijab(m,n,e,f)
!mc         AmpTerm = (-.5)*DubAmpl(j,n,f,b) + SingAmpl(j,f)*SingAmpl(n,b)
!mc         SumNF = SumNF + AmpTerm*ijab(m,n,e,f)
      End Do
   End Do
   Retval = SumNF + (-iajb(m,b,j,e)) + SumF + SumN
End Subroutine WBoth
!
!This program is designed to accept previously calculated Hartree-Fock two-electron integrals and apply the 
!    coupled-cluster post-SCF method. Currently, the program only incorporates contributions of Double excitations (CCD)
!This program will be carried out in a number of steps:
!    (1) Reading in the previously calculated integrals
!    (2) calculating an initial guess of the excitation amplitudes defined solely as the corresponding 
!          two-electron integrals
!    (3) iteratively solving for the excitation amplitudes and incorporating them into the final coupled-cluster Energy
!!This next block hoUses the actual program, which Calls the various Subroutines and calculates the amplitudes 
!    and Energy.
!
Subroutine CCD
   Use Global
   Implicit Integer (a-n)
   Double Precision, Allocatable, Dimension(:,:,:,:) :: UpdDubAmpl
   Character(300) :: FileName, Junk
   Double Precision ConvA, ConvE, Energy, PrevEnergy, Term1, Term2, Term3, Term4, Term5, Term6,&
&Term7, Term8, Term9, Term10, Retval, Ret
   Integer NO, NV, aaa, Counter
   Common NO, NV
!!First, the program Reads the previously calculated values based on occupied (i,j,k,l...m,n) orbitals and 
!    virtual (a,b,c,d...e,f) orbitals:
   Open(40, File='Input')
   Read(40,*)
   Read(40,*,IOStat=IError3) Junk, Junk, FileName
!EpsI  !default if it doesn't find the file is EpsI   
   Open(18, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(18,*) NO
     Allocate(EpsI(NO))
     Do i = 1, NO
        Read(18,*) EpsI(i)
     End Do
     IError2 = 0
   End If 
   Close(18)
   If (IError .ne. 0) Then
     Open(18, File='EpsI', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(18,*) NO
       Allocate(EpsI(NO))
       Do i = 1, NO
          Read(18,*) EpsI(i)
       End Do
     End If
     Close(18)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "EpsI not found"
     Return
   End If
!!!!!where we have Read the orbital energies (EpsI) of the Occupied orbitals
   Read(40,*,IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!EpsA  !default if it doesn't find the file is EpsA  
   Open(9, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(9,*) NV
     Allocate(EpsA(NV))
     Do i = 1, NV
        Read(9,*) EpsA(i)
     End Do
     IError2=0
   End If 
   Close(9)
   If (IError .ne. 0) Then
     Open(9, File='EpsA', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(9,*) NV
       Allocate(EpsA(NV))
       Do i = 1, NV
          Read(9,*) EpsA(i)
       End Do
     End If
     Close(9)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "EpsA not found"
     Return
   End If
!!!!!similarly, we have Read the energies of the Virtual orbitals
!!!!!Next, the two-electron integrals (of form <wx||yz>) are Read:
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iabcd  !default if it doesn't find the file is Iabcd 
   Open(11, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(11,*) NO, NV
     Allocate(abcd(NV,NV,NV,NV))
     Do d = 1, NV
        Do c = 1, NV
           Do b = 1, NV
              Do a = 1, NV
                 Read(11,*) abcd(a,b,c,d)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(11)
   If (IError .ne. 0) Then
     Open(11, File='Iabcd', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(11,*) NO, NV
       Allocate(abcd(NV,NV,NV,NV))
       Do d = 1, NV
          Do c = 1, NV
             Do b = 1, NV
                Do a = 1, NV
                   Read(11,*) abcd(a,b,c,d)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(11)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ab||cd> not found"
     Return
   End If
!!!!!here we have constructed the Global 4-Dimensional matrix of integrals whose indices run over 
!         only virtual orbitals (e.g. a,b,c,d)
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iiabc  !default if it doesn't find the file is Iiabc 
   Open(12, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(12,*) NO, NV
     Allocate(iabc(NO,NV,NV,NV))
     Do c = 1, NV
        Do b = 1, NV
           Do a = 1, NV
              Do i = 1, NO
                 Read(12,*) iabc(i,a,b,c)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(12)
   If (IError .ne. 0) Then
     Open(12, File='Iiabc', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(12,*) NO, NV
       Allocate(iabc(NO,NV,NV,NV))
       Do c = 1, NV
          Do b = 1, NV
             Do a = 1, NV
                Do i = 1, NO
                   Read(12,*) iabc(i,a,b,c)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(12)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ia||bc> not found"
     Return
   End If
!!!!!here we have constructed a similar matrix of integrals involving three virtual orbitals and one occupied orbital 
!
!!!!!similarly, the Global variable of two occupied, two virtual integrals in the form <ia||jb>:
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iiajb  !default if it doesn't find the file is Iiajb 
   Open(13, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(13,*) NO, NV
     Allocate(iajb(NO,NV,NO,NV))
     Do b = 1, NV
        Do j = 1, NO
           Do a = 1, NV
              Do i = 1, NO
                 Read(13,*) iajb(i,a,j,b)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(13)
   If (IError .ne. 0) Then
     Open(13, File='Iiajb', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(13,*) NO, NV
       Allocate(iajb(NO,NV,NO,NV))
       Do b = 1, NV
          Do j = 1, NO
             Do a = 1, NV
                Do i = 1, NO
                   Read(13,*) iajb(i,a,j,b)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(13)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ia||jb> not found"
     Return
   End If
!!!!!and the same for the form <ij||ab>:
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iijab  
   Open(14, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(14,*) NO, NV
     Allocate(ijab(NO,NO,NV,NV))
     Do b = 1, NV
        Do a = 1, NV
           Do j = 1, NO
              Do i = 1, NO
                 Read(14,*) ijab(i,j,a,b)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(14)
   If (IError .ne. 0) Then
     Open(14, File='Iijab', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(14,*) NO, NV
       Allocate(ijab(NO,NO,NV,NV))
       Do b = 1, NV
          Do a = 1, NV
             Do j = 1, NO
                Do i = 1, NO
                   Read(14,*) ijab(i,j,a,b)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(14)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ij||ab> not found"
     Return
   End If
!!!!!and for three occupied, one virtual (<ij||ka>):
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iijka  !default if it doesn't find the file is Iijka 
   Open(15, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(15,*) NO, NV
     Allocate(ijka(NO,NO,NO,NV))
     Do a = 1, NV
        Do k = 1, NO
           Do j = 1, NO
              Do i = 1, NO
                 Read(15,*) ijka(i,j,k,a)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(15)
   If (IError .ne. 0) Then
     Open(15, File='Iijka', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(15,*) NO, NV
       Allocate(ijka(NO,NO,NO,NV))
       Do a = 1, NV
          Do k = 1, NO
             Do j = 1, NO
                Do i = 1, NO
                   Read(15,*) ijka(i,j,k,a)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(15)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ij||ka> not found"
     Return
   End If
!!!!!and for integrals of only occupied orbitals:
   Read(40,*, IOStat=IError3) Junk, Junk, FileName
   IError2 = 400
!Iijkl  !default if it doesn't find the file is Iijkl 
   Open(16, File=FileName, Status='old', IOStat=IError)
   If (IError .eq. 0) Then
     Read(16,*) NO, NV
     Allocate(ijkl(NO,NO,NO,NO))
     Do l = 1, NO
        Do k = 1, NO
           Do j = 1, NO
              Do i = 1, NO
                 Read(16,*) ijkl(i,j,k,l)
              End Do
           End Do
        End Do
     End Do
     IError2 = 0
   End If 
   Close(16)
   If (IError .ne. 0) Then
     Open(16, File='Iijkl', Status='old', IOStat=IError2)
     If (IError2 .eq. 0) Then
       Read(16,*) NO, NV
       Allocate(ijkl(NO,NO,NO,NO))
       Do l = 1, NO
          Do k = 1, NO
             Do j = 1, NO
                Do i = 1, NO
                   Read(16,*) ijkl(i,j,k,l)
                End Do
             End Do
          End Do
       End Do
     End If
     Close(16)
   End If 
   If ((IError .ne. 0) .and. (IError2 .ne. 0)) Then
     Print*, "<ij||kl> not found"
     Return
   End If
   Close(40)
!!The next block of the main program constructs the amplitudes and ultimately deTermines the Energy. 
!!!!!First we Allocate memory to the 4-Dimensional matrices for amplitudes and the denominator Term 'Dijab' 
!         which (we see below) is the sum of occupied orbital energies minus the virtual orbitals they will be
!         excited to in the corresponding cluster excitation 
   Allocate(UpdDubAmpl(NO,NO,NV,NV))
   Allocate(DubAmpl(NO,NO,NV,NV))
   Allocate(Dijab(NO,NO,NV,NV))
   Allocate(AmplDiff(NO,NO,NV,NV))
!!!!!Thus we calculate the Dijab matrix:
   Do b = 1, NV
      Do a = 1, NV
         Do j = 1, NO
            Do i = 1, NO
               Dijab(i,j,a,b) = EpsI(i) + EpsI(j) - EpsA(a) - EpsA(b)
            End Do
         End Do
      End Do
   End Do
!!Here we perform our initial guess of the excitation amplitudes based solely on the two-electron integrals from
!    the SCF calculation: 
   Do b = 1, NV
      Do a = 1, NV
         Do j = 1, NO
            Do i = 1, NO
               DubAmpl(i,j,a,b) = ijab(i,j,a,b)/Dijab(i,j,a,b)
            End Do
         End Do
      End Do
   End Do
!!!!!We then perform an initial Energy calculation based on these amplitude guesses: 
   Do b = 1, NV
      Do a = 1, NV
         Do j = 1, NO
            Do i = 1, NO
               Energy = Energy + .25*DubAmpl(i,j,a,b)*ijab(i,j,a,b)
            End Do
         End Do
      End Do
   End Do
   print*, "MP2 Energy", Energy
!!The program now proceeds to perform the actual iterations.
!!!!!Each of the "Terms" is one of the parts of the Double-excitation cluster operator (T2) equation, as it is 
!         broken Down and explicitly defined before being recombined at the End.
   Energy = 0
   PrevEnergy = 0
   aaa = 0
   Counter = 0
!!!!!thus we have provided an initial guess of zero for each of our following variables, allowing for ease of summation
!!!!!The following block of code is the main structure of the iterative loop,including the test of 
!         and criteria for Convergence.
   Do
      If (Counter == 100) Then
         print*, "Convergence not achieved. Current energy after 100 cycles is ", Energy
         Exit
      End If
      If (aaa == 1) Then !Convergence
         Exit
      End If
!!!!!where we have initially deTermined that If the program Does not Converge within 100 iterations, it will cease.
      Do b = 1, NV
         Do a = 1, NV
            Do j = 1, NO
               Do i = 1, NO
                  Term1 = 0
                  Term2 = 0
                  Term3 = 0
                  Term4 = 0
                  Term5 = 0
                  Term6 = 0
                  Term7 = 0
                  Term8 = 0
                  Term9 = 0
                  Term10 = 0
!!!!!here we have constructed the iterative loops over all indices of occupied (i,j) and virtual (a,b) orbitals 
!         prior to actually defining the individual Terms
!!Terms 1 and 2 consist of two permutations of the virtual orbitals (ab) within a sum over a single virtual orbital
!    (e) of Terms of the "F" matrix [defined as a Subroutine] whose corresponding indices are virtual orbitals: 
                  Do e = 1, NV
                     Call FV(b,e,Retval)
                     Term1 = Term1 + DubAmpl(i,j,a,e) * Retval
                     Call FV(a,e,Retval)
                     Term2 = Term2 + DubAmpl(i,j,b,e) * Retval
                  End Do
!!!!!thus we have summed every (*,e) element of the virtual-orbital F matrix ['FVirt'], constructing the Terms 1 and 2
!         [F(b,e) and F(a,e) respectively] alongside the corresponding excitation amplitudes
!!Similarly, Terms 3 and 4 consist of the contributions where a single occupied orbital (m) dIffers in Terms containing
!    the occupied-orbital F matrix ['FOcc']; each Term is a permutation of the occupied orbitals (i,j): 
                  Do m = 1, NO
                     Call FO(m, j, Retval)
                     Term3 = Term3 + DubAmpl(i,m,a,b) * Retval
                     Call FO(m, i, Retval)
                     Term4 = Term4 + DubAmpl(j,m,a,b) * Retval
                  End Do
!!Terms 5 and 6 are simple sums over two occupied (m,n) and two virtual (e,f) orbitals respectively, and Call upon the
!    corresponding W matrices and excitation amplitudes: 
                  Do n = 1, NO
                     Do m = 1, NO
                        Call WO(m,n,i,j, Retval)
                        Term5 = Term5 + .5*DubAmpl(m,n,a,b)*Retval
                     End Do
                  End Do
                  Do f = 1, NV
                     Do e = 1, NV
                        Call WV(a,b,e,f, Retval)
                        Term6 = Term6 + .5*DubAmpl(i,j,e,f)*Retval
                     End Do
                  End Do
!!Terms 6-10 are the contributions of the mixed (two occupied, two virtual) W matrix ['WBoth'] summed over one dIffering
!    occupied and one dIffering virtual orbital; each Term is part of simultaneous permutations of (a,b) and (i,j):
                  Do e = 1, NV
                     Do m = 1, NO
                        Call WB(m,b,e,j, Retval)
                        Term7 = Term7 + DubAmpl(i,m,a,e)*Retval
                        Call WB(m,a,e,j, Retval)
                        Term8 = Term8 + DubAmpl(i,m,b,e)*Retval
                        Call WB(m,b,e,i, Retval)
                        Term9 = Term9 + DubAmpl(j,m,a,e)*Retval
                        Call WB(m,a,e,i, Retval)
                        Term10 = Term10 + DubAmpl(j,m,b,e)*Retval
                     End Do
                  End Do
!!Having defined all the Terms in the amplitude equation, we combine them all together in the 'updated' excitation
!    amplitude:
                  UpdDubAmpl(i,j,a,b) = (ijab(i,j,a,b)+Term1-Term2-Term3&
&+Term4+Term5+Term6+Term7-Term8-Term9+Term10)/Dijab(i,j,a,b)
               End Do
            End Do
         End Do
      End Do
      AmplDiff = UpdDubAmpl-DubAmpl      
      DubAmpl = UpdDubAmpl
!!!!!and thus the previous excitation amplitude is redefined as the updated value
!!Finally, the program recalculates the coupled-cluster Energy from the new amplitudes:
      Energy = 0
      Do b = 1, NV
         Do a = 1, NV
            Do j = 1, NO
               Do i = 1, NO
                  Energy = Energy + .25*DubAmpl(i,j,a,b)*ijab(i,j,a,b)
               End Do
            End Do
         End Do
      End Do
!!and checks for Convergence based on both the Energy dIfference (between current Energy and previous Energy), as well
!    as the dIfference in the root-mean squared (RMS) value of the excitation amplitudes:
      Call AmplConv(ConvA)
      ConvE = Energy - PrevEnergy
      ConvE = Sqrt(ConvE**2)
      If ((ConvA <= 1.0D-5) .and. (ConvE <= 1.0D-7)) then
         aaa = 1
      End If
!!!!!If Convergence has been achieved, the program signals the End of the iterative loop and prints the final Energy;
!         otherwise, it cycles back through the next iteration
      PrevEnergy = Energy
      Counter = Counter + 1
      Print*, "The energy of iteration ", Counter, "is", Energy
   End Do
   print*, "The final Energy is: ", Energy
End Subroutine CCD
!!The following blocks are each of the Subroutines/Terms incorporated into the overarching amplitude equation.
!
!!The F matrices correspond to their namesake and are defined by the nature of their indices - e.g. 'FV' 
!         is the F matrix defined only through the virtual orbitals of the excitation.
Subroutine FV(b, e, Retval)
   Use Global
   Implicit Integer (f-o)
   Integer, intent(in) :: b, e
   Double Precision, intent(out) :: Retval
   Double Precision Temp, TaoT
   Integer NO, NV
   Common NO, NV
   Temp = 0
   Do f = 1, NV
      Do n = 1, NO
         Do m = 1, NO
            Temp = Temp + DubAmpl(m,n,b,f)*ijab(m,n,e,f)
         End Do
      End Do
   End Do
   Retval = Temp * (-.5)
End Subroutine FV

Subroutine FO(m, j, Retval)
   Use Global
   Implicit None
   Integer, intent(in) :: j, m
   Double Precision, intent(out) :: Retval
   Double Precision Temp, TaoT
   Integer n, e, f, NO, NV
   Common NO, NV
   Temp = 0
   Do f = 1, NV
      Do e = 1, NV
         Do n = 1, NO
            Temp = Temp + DubAmpl(j,n,e,f)*ijab(m,n,e,f)
         End Do
      End Do
   End Do
   Retval = Temp * .5
End Subroutine FO
!
!!The W matrices similarly correspond to their namesake and are defined by the nature of their indices - e.g. 'WV' 
!         is the W matrix defined only through the virtual orbitals of the excitation.
Subroutine WO(m,n,i,j,Retval)
   Use Global
   Implicit None
   Integer, intent(in) :: m, n, i, j
   Double Precision, intent(out) :: Retval
   Double Precision Temp, Tao
   Integer e, f, NO, NV
   Common NO, NV
   Temp = 0
   Do f = 1, NV
      Do e = 1, NV
         Temp = Temp +.25*DubAmpl(i,j,e,f)*ijab(m,n,e,f)
      End Do
   End Do
   Retval = Temp + ijkl(m,n,i,j) 
End Subroutine WO
!
Subroutine WV(a,b,e,f,Retval)
   Use Global
   Implicit None
   Integer, intent(in) :: a, b, e, f
   Double Precision, intent(out) :: Retval
   Double Precision Temp, Tao
   Integer m, n, NO, NV
   Common NO, NV
   Temp = 0
   Do n = 1, NO
      Do m = 1, NO
         Temp = Temp +.25*DubAmpl(m,n,a,b)*ijab(m,n,e,f)
      End Do
   End Do
   Retval = Temp + abcd(a,b,e,f)
End Subroutine WV
!
Subroutine WB(m,b,e,j,Retval)
   Use Global
   Implicit None
   Integer, intent(in) :: m, b, e, j
   Double Precision, intent(out) :: Retval
   Double Precision Temp, Tao
   Integer n, f, NO, NV
   Common NO, NV
   Temp = 0
   Do f = 1, NV
      Do n = 1, NO
         Tao = (-.5)*DubAmpl(j,n,f,b)! + sampl(j,f)*sampl(n,b)
         Temp = Temp + Tao*ijab(m,n,e,f)
      End Do
   End Do
   Retval = Temp + (-iajb(m,b,j,e))
End Subroutine WB
!!The 'AmplConv' Subroutine is a simple method of computing the RMS of the excitation Amplitudes in order to obtain a
!    value ('Retval', passed into the main program as 'ConvA') with which to check Convergence.
Subroutine AmplConv(Retval)
   Use Global
   Implicit None
   Double Precision, Intent(out) :: Retval
   Integer a, b, i, j, NO, NV
   Double Precision Temp, m
   Common NO, NV
   m = NO * NO * NV * NV
   Temp = 0
   Do b = 1, NV
      Do a = 1, NV
         Do j = 1, NO
            Do i = 1, NO
               Temp = Temp + AmplDiff(i,j,a,b)**2
            End Do
         End Do
      End Do
   End Do
   Temp = Temp / m
   Retval = sqrt(Temp)
End Subroutine AmplConv
