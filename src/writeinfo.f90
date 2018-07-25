
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
      Character (10) :: dat, tim, acoord
      character*(77) :: string
      real(8) :: dumsum 

      acoord = "lattice"
      if (input%structure%cartesian) acoord = "cartesian"

      call printline(fnum,"=")
      write (string,'("EXCITING ", a, " started")') trim(versionname)
      call printtext(fnum,"=",string)
      if (len(trim(githash)) > 0) then
          write(string,'("version hash id: ",a)') githash
          call printtext(fnum,"=",string)
      end if
      call printtext(fnum,"=","")
#ifdef MPI
      Write (string,'("MPI version using ",i6," processor(s)")') procs
      call printtext(fnum,"=",string)
#ifndef MPI1   
      Write (string,'("|  using MPI-2 features")') 
      call printtext(fnum,"=",string)
#endif
#endif
 
      Call date_and_time (date=dat, time=tim)

      call printtext(fnum,"=","")
      Write (string,'("Date (DD-MM-YYYY) : ", A2, "-", A2, "-", A4)') &
     &  dat (7:8), dat (5:6), dat (1:4)
      call printtext(fnum,"=",string)
      Write (string,'("Time (hh:mm:ss)   : ", A2, ":", A2, ":", A2)') &
     &  tim (1:2), tim (3:4), tim (5:6)
      call printtext(fnum,"=",string)
      call printtext(fnum,"=","")
      Write (string,'("All units are atomic (Hartree, Bohr, etc.)")')
      call printtext(fnum,"=",string)
      call printline(fnum,"=")

      If (notelns .Gt. 0) Then
         Write (fnum,*)
         Write (fnum, '("Notes :")')
         Do i = 1, notelns
            Write (fnum, '(A)') notes (i)
         End Do
      End If

      Select Case (task)
      Case (0)
         call printbox(fnum,"*","Ground-state run starting from atomic densities")
      Case (1, 200)
         call printbox(fnum,"*","Ground-state run resuming from STATE.OUT")
      Case (2)
         call printbox(fnum,"*","Structural optimisation starting from atomic densities")
      Case (3)
         call printbox(fnum,"*","Structural optimisation run resuming from STATE.OUT")
      Case (5, 6)
         call printbox(fnum,"*","Ground-state Hartree-Fock run")
      Case (300)
         call printbox(fnum,"*","Reduced density-matrix functional theory run")
      Case Default
         Write (string,'("Error(writeinfo): task not defined : ", I8)') task
         call printtext(6,"#",string)
         Write (*,*)
         Stop
      End Select

      call printbox(fnum,"+","Starting initialization")

      Write (fnum,*)
      Write (fnum, '(" Lattice vectors (cartesian) :")')
      Write (fnum, '(3F18.10)') input%structure%crystal%basevect(1, 1), &
     & input%structure%crystal%basevect(2, 1), &
     & input%structure%crystal%basevect(3, 1)
      Write (fnum, '(3F18.10)') input%structure%crystal%basevect(1, 2), &
     & input%structure%crystal%basevect(2, 2), &
     & input%structure%crystal%basevect(3, 2)
      Write (fnum, '(3F18.10)') input%structure%crystal%basevect(1, 3), &
     & input%structure%crystal%basevect(2, 3), &
     & input%structure%crystal%basevect(3, 3)
      Write (fnum,*)
      Write (fnum, '(" Reciprocal lattice vectors (cartesian) :")')
      Write (fnum, '(3F18.10)') bvec (1, 1), bvec (2, 1), bvec (3, 1)
      Write (fnum, '(3F18.10)') bvec (1, 2), bvec (2, 2), bvec (3, 2)
      Write (fnum, '(3F18.10)') bvec (1, 3), bvec (2, 3), bvec (3, 3)
      Write (fnum,*)
      Write (fnum, '(" Unit cell volume", T45, ": ", F18.10)') omega
      Write (fnum, '(" Brillouin zone volume", T45, ": ", F18.10)') (twopi**3) / omega

      if (input%groundstate%outputlevelnumber>0) then
         If (input%structure%autormt) Then
            Write (fnum,*)
            Write (fnum, '(" Automatic determination of muffin-tin radii :")')
            Write (fnum, '("     parameters", T45, ": ", 2F18.10)') input%structure%rmtapm
         End If
         If (input%groundstate%frozencore) Then
            Write (fnum,*)
            Write (fnum, '(" A frozen-core calculation is performed")')
         End If
      end if

      dumsum = 0.
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            Do i = 1, 3
               dumsum = dumsum+abs(input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(i))
            End Do
         End Do
      End Do

      Do is = 1, nspecies
         Write (fnum,*)
         Write (fnum, '(" Species : ", I4, " (", A, ")")') is, trim (spsymb(is))
         Write (fnum, '("     parameters loaded from            ", T45, ":    ", A)') trim &
        & (input%structure%speciesarray(is)%species%speciesfile)
         Write (fnum, '("     name                              ", T45, ":    ", A)') trim (spname(is))
         if (input%groundstate%outputlevelnumber>0) then
            Write (fnum, '("     nuclear charge                    ", T45, ": ", F16.8)') spzn (is)
            Write (fnum, '("     electronic charge                 ", T45, ": ", F16.8)') spze (is)
            Write (fnum, '("     atomic mass                       ", T45, ": ", F16.8)') spmass (is)
            Write (fnum, '("     muffin-tin radius                 ", T45, ": ", F16.8)') rmt (is)
            Write (fnum, '("     # of radial points in muffin-tin  ", T45, ": ", I7)') nrmt (is) 
         end if
         Write (fnum,*)
         Write (fnum, '("     atomic positions (",A,") :")') trim(acoord)
         Do ia = 1, natoms (is)
            if (input%structure%cartesian) then  
                Write (fnum, '(T5,I4, " : ", 3F12.8)') ia, atposc(:,ia,is)
            else
                Write (fnum, '(T5,I4, " : ", 3F12.8)') ia, &
               &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            end if
         End Do

         if (dumsum .gt. input%structure%epslat) then
            Write (fnum,*)
            Write (fnum, '("     magnetic fields (cartesian) :")')
            Do ia = 1, natoms (is)
               Write (fnum, '(T5,I4, " : ", 3F12.8)') ia, &
              & input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
            End Do
         end if
      End Do

      Write (fnum,*)
      Write (fnum,    '(" Total number of atoms per unit cell   ", T45, ": ", I7)') natmtot
      Write (fnum,*)
      If (associated(input%groundstate%spin)) Then
         Write (fnum, '(" Spin treatment                        ", T45, ":    spin-polarised")')
      Else
         Write (fnum, '(" Spin treatment                        ", T45, ":    spin-unpolarised")')
      End If
      If (isspinorb()) Then
         Write (fnum, '(" Spin-orbit coupling is included")')
      End If
      If (associated(input%groundstate%spin)) Then
         Write (fnum, '(" Global magnetic field (cartesian) :",T42, 3F13.8)') input%groundstate%spin%bfieldc
         If (ncmag) Then
            Write (fnum, '(" non-collinear magnetisation")')
         Else
            Write (fnum, '(" collinear magnetisation in z-direction")')
         End If
      End If
      If (isspinspiral()) Then
         Write (fnum, '(" spin-spiral state assumed")')
         if (input%groundstate%outputlevelnumber>0) then
            Write (fnum, '("     q-vector (lattice)   : ", 3F16.8)') &
           & input%groundstate%spin%vqlss
            Write (fnum, '("     q-vector (cartesian) : ", 3F16.8)') vqcss
            Write (fnum, '("     q-vector length      : ", F16.8)') &
           & Sqrt(vqcss(1)**2+vqcss(2)**2+vqcss(3)**2)
         end if
      End If
      If (getfixspinnumber() .Ne. 0) Then
         Write (fnum, '(" fixed spin moment (FSM) calculation")')
      End If
      If ((getfixspinnumber() .Eq. 1) .Or. (getfixspinnumber() .Eq. 3)) Then
         Write (fnum, '("     fixed total spin moment (cartesian) : ", 3F16.8 )') &
        &  input%groundstate%spin%momfix
      End If
      If ((getfixspinnumber() .Eq. 2) .Or. (getfixspinnumber() .Eq. 3)) Then
         Write (fnum, '("     fixing local MT moments (cartesian) : ")')
         Do is = 1, nspecies
            Write (fnum, '("  species : ", I4, " (", A, ")")') is, trim (spsymb(is))
            Do ia = 1, natoms (is)
               Write (fnum, '("	", I4, 3G18.10)') &
              &  ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%mommtfix(:)
            End Do
         End Do
      End If

      Write (fnum,*)
      Write (fnum, '(" Number of Bravais lattice symmetries  ", T45, ": ", I7)') nsymlat
      Write (fnum, '(" Number of crystal symmetries          ", T45, ": ", I7)') nsymcrys

      If (input%groundstate%autokpt) Then
         Write (fnum,*)
         Write (fnum, '(" Sphere radius used for k-point grid density ", T45, ": ", F16.8)')&
        & input%groundstate%radkpt
      End If
      Write (fnum,*)
      Write (fnum, '(" k-point grid                          ", T45, ": ", I7,2I5)') input%groundstate%ngridk
      dumsum = 0.
      Do i = 1, 3
         dumsum = dumsum+abs(input%groundstate%vkloff(i))
      End Do
      if (dumsum .gt. input%structure%epslat) then
         Write (fnum, '(" k-point offset :     ", 3F10.6)') &
        & input%groundstate%vkloff
      end if
      Write (fnum, '(" Total number of k-points              ", T45, ": ", I7)') nkpt
      If (input%groundstate%reducek) Then
         Write (fnum, '(" k-point set is reduced with crystal symmetries")')
      Else
         Write (fnum, '(" k-point set is not reduced")')
      End If

      Write (fnum,*)
      Write (fnum, '(" R^MT_min * |G+k|_max (rgkmax)",T45,": ", F16.8)') input%groundstate%rgkmax
      if (input%groundstate%outputlevelnumber>0) then
         If ((input%groundstate%isgkmax .Ge. 1) .And. (input%groundstate%isgkmax .Le. nspecies)) Then
            Write (fnum, '(" Species with R^MT_min",T45,": ", I7, " (", A, ")")') &
           & input%groundstate%isgkmax, trim (spsymb(input%groundstate%isgkmax))
         End If
         Write (fnum, '(" Maximum |G+k| for APW functions",T45,": ", F16.8)') gkmax
         Write (fnum, '(" Maximum |G| for potential and density",T45,": ", F16.8)') &
        & input%groundstate%gmaxvr
      end if
      if (input%groundstate%outputlevelnumber>1) &
     & Write (fnum, '(" Polynomial order for pseudochg. density",T45,": ", I7)') &
     & input%groundstate%npsden

      Write (fnum,*)
      Write (fnum, '(" G-vector grid sizes                   ", T45, ":  ", 3I6)') ngrid (1), ngrid &
     & (2), ngrid (3)
      Write (fnum, '(" Total number of G-vectors             ", T45, ":", I8)') ngvec

      if (input%groundstate%outputlevelnumber>0) Then
         Write (fnum,*)
         Write (fnum, '(" Maximum angular momentum used for")')
         Write (fnum, '("     APW functions                     ", T45, ": ", I7)') &
        & input%groundstate%lmaxapw
         Write (fnum, '("     computing H and O matrix elements ", T45, ": ", I7)') &
        & input%groundstate%lmaxmat
         Write (fnum, '("     potential and density             ", T45, ": ", I7)') &
        & input%groundstate%lmaxvr
         Write (fnum, '("     inner part of muffin-tin          ", T45, ": ", I7)') &
        & input%groundstate%lmaxinr
      end if

      Write (fnum,*)
      Write (fnum, '(" Total nuclear charge                  ", T45, ": ", F16.8)') chgzn
      Write (fnum, '(" Total electronic charge               ", T45, ": ", F16.8)') chgtot
      if (input%groundstate%outputlevelnumber>0) Then
         Write (fnum, '(" Total core charge                     ", T45, ": ", F16.8)') chgcr
         Write (fnum, '(" Total valence charge                  ", T45, ": ", F16.8)') chgval
         if (abs(input%groundstate%chgexs)>input%structure%epslat) &
        & Write (fnum, '(" Total excess charge                   ", T45, ": ", F16.8)') &
        & input%groundstate%chgexs
      end if

      if (input%groundstate%outputlevelnumber>1) Then
         Write (fnum,*)
         Write (fnum, '(" Effective Wigner radius, r_s          ", T45, ": ", F16.8)') rwigner
      end if

      Write (fnum,*)
      Write (fnum, '(" Number of empty states                ", T45, ": ", I7)') &
     & input%groundstate%nempty
      Write (fnum, '(" Total number of valence states        ", T45, ": ", I7)') nstsv
      if (input%groundstate%outputlevelnumber>0) Then
         Write (fnum,*)
         Write (fnum, '(" Maximum Hamiltonian size      ", T45, ": ", I7)') nmatmax
         Write (fnum, '(" Maximum number of plane-waves ", T45, ": ", I7)') ngkmax
         Write (fnum, '(" Total number of local-orbitals", T45, ": ", I7)') nlotot
      end if

      Write (fnum,*)
      If ((task .Eq. 5) .Or. (task .Eq. 6)) &
     & write (fnum, '(" Hartree-Fock calculation using Kohn-Sham states")')
      If (input%groundstate%xctypenumber .Lt. 0) Then
         Write (fnum,*)
         Write (fnum, '(" Optimised effective potential (OEP) and exact exchange (EXX)")')
         Write (fnum, '("     Phys. Rev. B 53, 7024 (1996)")')
         Write (fnum, '("     Correlation type ", T45, ": ", I7)') Abs(input%groundstate%xctypenumber)
         Write (fnum, '("     ", A)') trim (xcdescr)
      Else
         Write (fnum, '(" Exchange-correlation type ", T45, ": ", I7)') input%groundstate%xctypenumber
         Write (fnum, '("     ", A)') trim (xcdescr)
      End If
      If (associated(input%groundstate%Hybrid )) Then
         write(fnum,  '("     Hybrid functional ")')
         write(fnum,  '("     Exchange type                         ", T45, ": ", A8)') trim (input%groundstate%Hybrid%exchangetype)
         write(fnum,  '("     Mixing coefficient for exact exchange ", T45, ": ", F16.8)') input%groundstate%Hybrid%excoeff
      Else If (xcgrad .Eq. 1)  Then 
         write (fnum, '("     Generalised gradient approximation (GGA)")')
      End If
      If (ldapu .Ne. 0) Then
         Write (fnum,*)
         Write (fnum, '(" LDA+U calculation")')
         If (ldapu .Eq. 1) Then
            Write (fnum, '("     fully localised limit (FLL)")')
         Else If (ldapu .Eq. 2) Then
            Write (fnum, '("     around mean field (AFM)")')
         Else If (ldapu .Eq. 3) Then
            Write (fnum, '("     interpolation between FLL and AFM")')
         Else
            Write (*,*)
            Write (*, '(" Error(writeinfo): ldapu not defined : ", I8)') ldapu
            Write (*,*)
            Stop
         End If
         Write (fnum, '("     see PRB 67, 153106 (2003) and PRB 52, R5467 (1995)")')
         if (input%groundstate%outputlevelnumber>0) Then
            Do is = 1, nspecies
               If (llu(is) .Ge. 0) Then
                  Write (fnum, '("     species : ", I4, " (", A, ")", ", l = "&
                 &, I2, ", U = ", F12.8, ", J = ", F12.8)') is, trim &
                 & (spsymb(is)), llu (is), ujlu (1, is), ujlu (2, is)
               End If
            End Do
         end if
      End If
      If (task .Eq. 300) Then
         Write (fnum,*)
         Write (fnum, '(" RDMFT calculation")')
         Write (fnum, '("     see arXiv:0801.3787v1 [cond-mat.mtrl-sci]")')
         Write (fnum, '("     RDMFT exchange-correlation type ", T45, ": ", I7)') &
        & input%groundstate%RDMFT%rdmxctype
         If (input%groundstate%RDMFT%rdmxctype .Eq. 1) Then
            Write (fnum, '("     Hartree-Fock functional")')
         Else If (input%groundstate%RDMFT%rdmxctype .Eq. 2) Then
            Write (fnum, '("     SDLG functional, exponent : ", F16.8)') &
           & input%groundstate%RDMFT%rdmalpha
         End If
      End If

      Write (fnum,*)
#ifdef TETRA
      tetocc=.false.
      If (associated(input%xs)) Then
         If (associated(input%xs%tetra)) Then
            tetocc=input%xs%tetra%tetraocc
         end if
      end if
      If ( .Not. tetocc) Then
#endif
         Write (fnum, '(" Smearing scheme", T45, ":    ", A)') trim (sdescr)
         if (input%groundstate%stypenumber.ge.0) then
           Write (fnum, '(" Smearing width", T45, ": ", F16.8)') &
          & input%groundstate%swidth
         end if
#ifdef TETRA
      Else
         Write (fnum, '(" Using linear tetrahedron method ")')
      End If
#endif
  
      if (input%groundstate%outputlevelnumber>0) Then
         Write (fnum,*)
         select case(input%groundstate%mixernumber)
           case(1)
             write(fnum,'(" Using adaptive step size linear potential mixing")')
           case(2)
             write(fnum,'(" Using multisecant Broyden potential mixing")')
           case(3)
             write(fnum,'(" Using Pulay potential mixing")')
         end select
      end if 

      call printbox(fnum,"+","Ending initialization")

      Call flushifc (fnum)
      Return
End Subroutine
!EOC
