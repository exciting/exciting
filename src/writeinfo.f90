
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeinfo
! !INTERFACE:
!
!
Subroutine writeinfo (fnum)
! !USES:
      Use modinput
      Use modmain
      use modmpi, only: procs
#ifdef TETRA
      Use modtetra
#endif
! !INPUT/OUTPUT PARAMETERS:
!   fnum : unit specifier for INFO.OUT file (in,integer)
! !DESCRIPTION:
!   Outputs basic information about the run to the file {\tt INFO.OUT}. Does not
!   close the file afterwards.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer :: fnum
! local variables
      Integer :: i, is, ia
#ifdef TETRA
      logical :: tetocc
#endif
      Character (10) :: dat, tim
      Write (fnum, '("+-----------------------------------------------------------+")')
      Write (fnum, '("| EXCITING helium    (",I2.2,".",I2.2,".",I2.2,") started                     |")') version
      Write (fnum, '("| version hash id: ",a," |")') githash
#ifdef MPI
      Write (fnum, '("| MPI version using ",i6," processor(s)                     |")') procs
#ifndef MPI1
      Write (fnum, '("|  using MPI-2 features                                     |")')
#endif
#endif
      Write (fnum, '("+-----------------------------------------------------------+")')
      If (notelns .Gt. 0) Then
         Write (fnum,*)
         Write (fnum, '("Notes :")')
         Do i = 1, notelns
            Write (fnum, '(A)') notes (i)
         End Do
      End If
      Call date_and_time (date=dat, time=tim)
      Write (fnum,*)
      Write (fnum, '("Date (YYYY-MM-DD) : ", A4, "-", A2, "-", A2)') &
     & dat (1:4), dat (5:6), dat (7:8)
      Write (fnum, '("Time (hh:mm:ss)   : ", A2, ":", A2, ":", A2)') tim &
     & (1:2), tim (3:4), tim (5:6)
      Write (fnum,*)
      Write (fnum, '("All units are atomic (Hartree, Bohr, etc.)")')
      Select Case (task)
      Case (0)
         Write (fnum,*)
         Write (fnum, '("+---------------------------------------------&
        &----+")')
         Write (fnum, '("| Ground-state run starting from atomic densit&
        &ies |")')
         Write (fnum, '("+---------------------------------------------&
        &----+")')
      Case (1, 200)
         Write (fnum,*)
         Write (fnum, '("+------------------------------------------+")&
        &')
         Write (fnum, '("| Ground-state run resuming from STATE.OUT |")&
        &')
         Write (fnum, '("+------------------------------------------+")&
        &')
      Case (2)
         Write (fnum,*)
         Write (fnum, '("+---------------------------------------------&
        &-----------+")')
         Write (fnum, '("| Structural optimisation starting from atomic&
        & densities |")')
         Write (fnum, '("+---------------------------------------------&
        &-----------+")')
      Case (3)
         Write (fnum,*)
         Write (fnum, '("+-----------------------------------------------------+")')
         Write (fnum, '("| Structural optimisation run resuming from STATE.OUT |")')
         Write (fnum, '("+-----------------------------------------------------+")')
      Case (5, 6)
         Write (fnum,*)
         Write (fnum, '("+-------------------------------+")')
         Write (fnum, '("| Ground-state Hartree-Fock run |")')
         Write (fnum, '("+-------------------------------+")')
      Case (300)
         Write (fnum,*)
         Write (fnum, '("+---------------------------------------------&
        &-+")')
         Write (fnum, '("| Reduced density matrix functional theory run&
        & |")')
         Write (fnum, '("+---------------------------------------------&
        &-+")')
      Case Default
         Write (*,*)
         Write (*, '("Error(writeinfo): task not defined : ", I8)') &
        & task
         Write (*,*)
         Stop
      End Select
      Write (fnum,*)
      Write (fnum, '("Lattice vectors :")')
      Write (fnum, '(3G18.10)') input%structure%crystal%basevect(1, 1), &
     & input%structure%crystal%basevect(2, 1), &
     & input%structure%crystal%basevect(3, 1)
      Write (fnum, '(3G18.10)') input%structure%crystal%basevect(1, 2), &
     & input%structure%crystal%basevect(2, 2), &
     & input%structure%crystal%basevect(3, 2)
      Write (fnum, '(3G18.10)') input%structure%crystal%basevect(1, 3), &
     & input%structure%crystal%basevect(2, 3), &
     & input%structure%crystal%basevect(3, 3)
      Write (fnum,*)
      Write (fnum, '("Reciprocal lattice vectors :")')
      Write (fnum, '(3G18.10)') bvec (1, 1), bvec (2, 1), bvec (3, 1)
      Write (fnum, '(3G18.10)') bvec (1, 2), bvec (2, 2), bvec (3, 2)
      Write (fnum, '(3G18.10)') bvec (1, 3), bvec (2, 3), bvec (3, 3)
      Write (fnum,*)
      Write (fnum, '("Unit cell volume      : ", G18.10)') omega
      Write (fnum, '("Brillouin zone volume : ", G18.10)') (twopi**3) / &
     & omega
      If (input%structure%autormt) Then
         Write (fnum,*)
         Write (fnum, '("Automatic determination of muffin-tin radii")')
         Write (fnum, '(" parameters : ", 2G18.10)') &
        & input%groundstate%rmtapm
      End If
      If (input%groundstate%frozencore) Then
         Write (fnum,*)
         Write (fnum, '("A frozen-core calculation is performed")')
      End If
      Do is = 1, nspecies
         Write (fnum,*)
         Write (fnum, '("Species : ", I4, " (", A, ")")') is, trim &
        & (spsymb(is))
         Write (fnum, '(" parameters loaded from : ", A)') trim &
        & (input%structure%speciesarray(is)%species%speciesfile)
         Write (fnum, '(" name : ", A)') trim (spname(is))
         Write (fnum, '(" nuclear charge    : ", G18.10)') spzn (is)
         Write (fnum, '(" electronic charge : ", G18.10)') spze (is)
         Write (fnum, '(" atomic mass : ", G18.10)') spmass (is)
         Write (fnum, '(" muffin-tin radius : ", G18.10)') rmt (is)
         Write (fnum, '(" number of radial points in muffin-tin : ", I6&
        &)') nrmt (is)
         Write (fnum, '(" atomic positions (lattice), magnetic fields (&
        &Cartesian) :")')
         Do ia = 1, natoms (is)
            Write (fnum, '(I4, " : ", 3F12.8, "	", 3F12.8)') ia, &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
         End Do
      End Do
      Write (fnum,*)
      Write (fnum, '("Total number of atoms per unit cell : ", I4)') &
     & natmtot
      Write (fnum,*)
      Write (fnum, '("Spin treatment :")')
      If (associated(input%groundstate%spin)) Then
         Write (fnum, '(" spin-polarised")')
      Else
         Write (fnum, '(" spin-unpolarised")')
      End If
      If (isspinorb()) Then
         Write (fnum, '(" spin-orbit coupling")')
      End If
      If (associated(input%groundstate%spin)) Then
         Write (fnum, '(" global magnetic field (Cartesian) : ", 3G18.1&
        &0)') input%groundstate%spin%bfieldc
         If (ncmag) Then
            Write (fnum, '(" non-collinear magnetisation")')
         Else
            Write (fnum, '(" collinear magnetisation in z-direction")')
         End If
      End If
      If (isspinspiral()) Then
         Write (fnum, '(" spin-spiral state assumed")')
         Write (fnum, '("  q-vector (lattice)	: ", 3G18.10)') &
        & input%groundstate%spin%vqlss
         Write (fnum, '("  q-vector (Cartesian) : ", 3G18.10)') vqcss
         Write (fnum, '("  q-vector length	: ", G18.10)') Sqrt &
        & (vqcss(1)**2+vqcss(2)**2+vqcss(3)**2)
      End If
      If (getfixspinnumber() .Ne. 0) Then
         Write (fnum, '(" fixed spin moment (FSM) calculation")')
      End If
      If ((getfixspinnumber() .Eq. 1) .Or. (getfixspinnumber() .Eq. 3)) &
     & Then
         Write (fnum, '("  fixing total moment to (Cartesian) :")')
         Write (fnum, '("  ", 3G18.10)') input%groundstate%spin%momfix
      End If
      If ((getfixspinnumber() .Eq. 2) .Or. (getfixspinnumber() .Eq. 3)) &
     & Then
         Write (fnum, '("  fixing local muffin-tin moments to (Cartesia&
        &n) :")')
         Do is = 1, nspecies
            Write (fnum, '("  species : ", I4, " (", A, ")")') is, trim &
           & (spsymb(is))
            Do ia = 1, natoms (is)
               Write (fnum, '("	", I4, 3G18.10)') ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%mommtfix(:)
            End Do
         End Do
      End If
      Write (fnum,*)
      Write (fnum, '("Number of Bravais lattice symmetries : ", I4)') &
     & nsymlat
      Write (fnum, '("Number of crystal symmetries	    : ", I4)') &
     & nsymcrys
      Write (fnum,*)
      If (input%groundstate%autokpt) Then
         Write (fnum, '("radius of sphere used to determine k-point gri&
        &d density : ", G18.10)') input%groundstate%radkpt
      End If
      Write (fnum, '("k-point grid : ", 3I6)') input%groundstate%ngridk
      Write (fnum, '("k-point offset : ", 3G18.10)') &
     & input%groundstate%vkloff
      If (input%groundstate%reducek) Then
         Write (fnum, '("k-point set is reduced with crystal symmetries&
        &")')
      Else
         Write (fnum, '("k-point set is not reduced")')
      End If
      Write (fnum, '("Total number of k-points : ", I8)') nkpt
      Write (fnum,*)
      Write (fnum, '("Smallest muffin-tin radius times maximum |G+k| : &
     &", G18.10)') input%groundstate%rgkmax
      If ((input%groundstate%isgkmax .Ge. 1) .And. &
     & (input%groundstate%isgkmax .Le. nspecies)) Then
         Write (fnum, '("Species with smallest (or selected) muffin-tin&
        & radius : ", I4, " (", A, ")")') input%groundstate%isgkmax, &
        & trim (spsymb(input%groundstate%isgkmax))
      End If
      Write (fnum, '("Maximum |G+k| for APW functions	     : ", G18.10)&
     &') gkmax
      Write (fnum, '("Maximum |G| for potential and density : ", G18.10&
     &)') input%groundstate%gmaxvr
      Write (fnum, '("Polynomial order for pseudocharge density : ", I4&
     &)') input%groundstate%npsden
      Write (fnum,*)
      Write (fnum, '("G-vector grid sizes : ", 3I6)') ngrid (1), ngrid &
     & (2), ngrid (3)
      Write (fnum, '("Total number of G-vectors : ", I8)') ngvec
      Write (fnum,*)
      Write (fnum, '("Maximum angular momentum used for")')
      Write (fnum, '(" APW functions                     : ", I4)') &
     & input%groundstate%lmaxapw
      Write (fnum, '(" computing H and O matrix elements : ", I4)') &
     & input%groundstate%lmaxmat
      Write (fnum, '(" potential and density             : ", I4)') &
     & input%groundstate%lmaxvr
      Write (fnum, '(" inner part of muffin-tin          : ", I4)') &
     & input%groundstate%lmaxinr
      Write (fnum,*)
      Write (fnum, '("Total nuclear charge    : ", G18.10)') chgzn
      Write (fnum, '("Total core charge       : ", G18.10)') chgcr
      Write (fnum, '("Total valence charge    : ", G18.10)') chgval
      Write (fnum, '("Total excess charge     : ", G18.10)') &
     & input%groundstate%chgexs
      Write (fnum, '("Total electronic charge : ", G18.10)') chgtot
      Write (fnum,*)
      Write (fnum, '("Effective Wigner radius, r_s : ", G18.10)') &
     & rwigner
      Write (fnum,*)
      Write (fnum, '("Number of empty states	      : ", I4)') &
     & input%groundstate%nempty
      Write (fnum, '("Total number of valence states : ", I4)') nstsv
      Write (fnum,*)
      Write (fnum, '("Total number of local-orbitals : ", I4)') nlotot
      Write (fnum,*)
      If ((task .Eq. 5) .Or. (task .Eq. 6)) write (fnum, '("Hartree-Foc&
     &k calculation using Kohn-Sham states")')
      If (input%groundstate%xctypenumber .Lt. 0) Then
         Write (fnum, '("Optimised effective potential (OEP) and exact &
        &exchange (EXX)")')
         Write (fnum, '(" Phys. Rev. B 53, 7024 (1996)")')
         Write (fnum, '("Correlation type : ", I4)') Abs &
        & (input%groundstate%xctypenumber)
         Write (fnum, '(" ", A)') trim (xcdescr)
      Else
         Write (fnum, '("Exchange-correlation type : ", I4)') &
        & input%groundstate%xctypenumber
         Write (fnum, '(" ", A)') trim (xcdescr)
      End If
      If (xcgrad .Eq. 1) write (fnum, '(" Generalised gradient approxim&
     &ation (GGA)")')
      If (ldapu .Ne. 0) Then
         Write (fnum,*)
         Write (fnum, '("LDA+U calculation")')
         If (ldapu .Eq. 1) Then
            Write (fnum, '(" fully localised limit (FLL)")')
         Else If (ldapu .Eq. 2) Then
            Write (fnum, '(" around mean field (AFM)")')
         Else If (ldapu .Eq. 3) Then
            Write (fnum, '(" interpolation between FLL and AFM")')
         Else
            Write (*,*)
            Write (*, '("Error(writeinfo): ldapu not defined : ", I8)') &
           & ldapu
            Write (*,*)
            Stop
         End If
         Write (fnum, '(" see PRB 67, 153106 (2003) and PRB 52, R5467 (&
        &1995)")')
         Do is = 1, nspecies
            If (llu(is) .Ge. 0) Then
               Write (fnum, '(" species : ", I4, " (", A, ")", ", l = "&
              &, I2, ", U = ", F12.8, ", J = ", F12.8)') is, trim &
              & (spsymb(is)), llu (is), ujlu (1, is), ujlu (2, is)
            End If
         End Do
      End If
      If (task .Eq. 300) Then
         Write (fnum,*)
         Write (fnum, '("RDMFT calculation")')
         Write (fnum, '(" see arXiv:0801.3787v1 [cond-mat.mtrl-sci]")')
         Write (fnum, '(" RDMFT exchange-correlation type : ", I4)') &
        & input%groundstate%RDMFT%rdmxctype
         If (input%groundstate%RDMFT%rdmxctype .Eq. 1) Then
            Write (fnum, '("  Hartree-Fock functional")')
         Else If (input%groundstate%RDMFT%rdmxctype .Eq. 2) Then
            Write (fnum, '("  SDLG functional, exponent : ", G18.10)') &
           & input%groundstate%RDMFT%rdmalpha
         End If
      End If
      Write (fnum,*)
      Write (fnum, '("Smearing scheme :")')
#ifdef TETRA
      tetocc=.false.
      If (associated(input%xs)) Then
         If (associated(input%xs%tetra)) Then
            tetocc=input%xs%tetra%tetraocc
         end if
      end if
      If ( .Not. tetocc) Then
#endif
         Write (fnum, '(" ", A)') trim (sdescr)
         Write (fnum, '("Smearing width : ", G18.10)') &
          & input%groundstate%swidth
#ifdef TETRA
      Else
         Write (fnum, '(" ", A)') 'No smearing - using the linear&
         & tetrahedron method'
         Write (fnum, '(" ", A)') 'for occupation numbers and Fer&
         &mi energy'
      End If
#endif
      Write (fnum,*)
      Write (fnum, '("Radial integration step length : ", I4)') &
     & input%groundstate%lradstep
      Call flushifc (fnum)
      Return
End Subroutine
!EOC
