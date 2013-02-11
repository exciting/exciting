! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gndstate
! !INTERFACE:
!
!
Subroutine gndstate
! !USES:
    Use modinput
    Use modmain
    Use modmpi
    Use scl_xml_out_Module
!
! !DESCRIPTION:
!   Computes the self-consistent Kohn-Sham ground-state. General information is
!   written to the file {\tt INFO.OUT}. First- and second-variational
!   eigenvalues, eigenvectors and occupancies are written to the unformatted
!   files {\tt EVALFV.OUT}, {\tt EVALSV.OUT}, {\tt EVECFV.OUT}, {\tt EVECSV.OUT}
!   and {\tt OCCSV.OUT}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
    Implicit None

! time measurements
    Real(8) :: timetot, ts0, ts1, tsg0, tsg1, tin1, tin0

!! TIME - Initialisation segment
    Call timesec (tsg0)
    Call timesec (ts0)
     
! initialise global variables
    Call timesec (tin0)
    Call init0
    Call timesec (tin1)
    time_init0=tin1-tin0
    Call timesec (tin0)
    Call init1
    Call timesec (tin1)
    time_init1=tin1-tin0
    
! initialise OEP variables if required
    if ((input%groundstate%xctypenumber .Lt. 0) .Or. &
   &    (xctype(2) .Ge. 400) .Or. (xctype(1) .Ge. 400)) then
      call init2
    end if
    If (rank .Eq. 0) Then
! write the real and reciprocal lattice vectors to file
      Call writelat
! write interatomic distances to file
      Call writeiad (.False.)
! write symmetry matrices to file
      Call writesym
#ifdef XS
! write realtion to inverse symmetries
      Call writesymi
! write advanced information on symmetry group
      Call writesym2
! write out symmetrization matrix for rank 2 tensors
      Call writesymt2
#endif
! output the k-point set to file
      Call writekpts
! write lattice vectors and atomic positions to file
      Call writegeom (.False.)
      Call writegeometryxml (.False.)
! open INFO.OUT file
      Open (60, File='INFO'//trim(filext), Action='WRITE', Form='FORMATTED')
! open TOTENERGY.OUT
      Open (61, File='TOTENERGY'//trim(filext), Action='WRITE', Form='FORMATTED')
! open FERMIDOS.OUT
      Open (62, File='FERMIDOS'//trim(filext), Action='WRITE', Form='FORMATTED')
! open MOMENT.OUT if required
      If (associated(input%groundstate%spin)) open (63, file='MOMENT'//trim(filext), action='WRITE', form='FORMATTED')
! open FORCEMAX.OUT if required
      If (input%groundstate%tforce) open (64, file='FORCEMAX'//trim(filext), action='WRITE', form='FORMATTED')
! open RMSDVEFF.OUT
      Open (65, File='RMSDVEFF'//trim(filext), Action='WRITE', Form='FORMATTED')
! open DTOTENERGY.OUT
      open(66,file='DTOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open DFORCEMAX.OUT
      if (input%groundstate%tforce) open(67,file='DFORCEMAX'//trim(filext),action='WRITE',form='FORMATTED')
! open CHGDIST.OUT
      open(68,file='CHGDIST'//trim(filext),action='WRITE',form='FORMATTED')
! open PCHARGE.OUT
      open(69,file='PCHARGE'//trim(filext),action='WRITE',form='FORMATTED')
! write out general information to INFO.OUT
      Call writeinfo (60)
      write (60,*)
    End If
! delete any existing eigenvector files
    If ((rank .Eq. 0) .And. ((task .Eq. 0) .Or. (task .Eq. 2))) Call delevec
    Call timesec (ts1)
    timeinit = timeinit+ts1-ts0
!! TIME - End of initialisation segment

!------------------------------------!
!   SCF cycle
!------------------------------------!
    call scf_cycle
 
!------------------------------------!
!   perform structure optimization
!------------------------------------!
    If (( .Not. tstop) .And. ((task .Eq. 2) .Or. (task .Eq. 3))) Then

        If (rank .Eq. 0) Then
            Write (60,*)
            Write (60, '("Maximum force magnitude (target) : ", G18.10, &
           &  " (", G18.10, ")")') forcemax, input%structureoptimization%epsforce
            Call flushifc (60)
        End If

        ! check force convergence
        If (forcemax .Le. input%structureoptimization%epsforce) Then
            If (rank .Eq. 0) Then
                Write (60,*)
                Write (60, '("Force convergence target achieved")')
            End If
        Else
            call structureoptimization()
        End If

    End If

!! TIME - Fifth IO segment
    Call timesec(ts0)
!-------------------------------------!
!   output timing information
!-------------------------------------!
    If (rank .Eq. 0) Then
        Write (60,*)
        Write (60, '("Timings (CPU seconds) :")')
        Write (60, '(" initialisation", T40, ": ", F12.2)') timeinit
        Write (60, '(" Hamiltonian and overlap matrix set up", T40,": ", F12.2)') timemat
        Write (60, '(" first-variational secular equation", T40, ": ", F12.2)') timefv
        If (associated(input%groundstate%spin)) Then
          Write (60, '(" second-variational calculation", T40, ": ", F12.2)') timesv
        End If
        Write (60, '(" charge density calculation", T40, ": ", F12.2)') timerho
        Write (60, '(" potential calculation", T40, ": ", F12.2)') timepot
        If (input%groundstate%tforce) Then
          Write (60, '(" force calculation", T40, ": ", F12.2)') timefor
        End If
        timetot = timeinit + timemat + timefv + timesv + timerho + timepot + timefor
        Write (60, '(" total", T40, ": ", F12.2)') timetot
        Write (60,*)
        Write (60, '("+-------------------------+")')
        Write (60, '("| EXCITING Lithium stopped |")')
        Write (60, '("+-------------------------+")')
! close the TOTENERGY.OUT file
        Close (61)
! close the FERMIDOS.OUT file
        Close (62)
! close the MOMENT.OUT file
        If (associated(input%groundstate%spin)) close (63)
! close the FORCEMAX.OUT file
        If (input%groundstate%tforce) close (64)
! close the RMSDVEFF.OUT file
        Close (65)
! close the DTOTENERGY.OUT file
        close(66)
! close the DFORCEMAX.OUT file
        if (input%groundstate%tforce) close(67)
! close the CHGDIST.OUT file
        close(68)
! close the PCHARGE.OUT file
        close(69)
        Call structure_xmlout ()
        Call scl_xml_setGndstateStatus ("finished")
        Call scl_xml_out_write ()
    End If
    Call timesec(ts1)
    timeio=ts1-ts0+timeio
!! TIME - End of fifth IO segment

    If (rank .Eq. 0) then
        Write (60,*)
        Write (60, '("Timings (seconds) :")')
        Write (60, '(" initialisation", T40, ": ", F12.2)') timeinit
        Write (60, '("            - init0", T40,": ", F13.2)') time_init0
        Write (60, '("            - init1", T40,": ", F13.2)') time_init1
        Write (60, '("            - rhoinit", T40,": ", F13.2)') time_density_init
        Write (60, '("            - potential initialisation", T40,": ", F13.2)') time_pot_init
        Write (60, '("            - others", T40,": ", F13.2)') timeinit-time_init0-time_init1
        Write (60, '(" Hamiltonian and overlap matrix set up", T40,&
       &  ": ", F12.2)') timemat
        Write (60, '("            - hmlaan", T40,": ", F13.2)') time_hmlaan
        Write (60, '("            - hmlalon", T40,": ", F13.2)') time_hmlalon
        Write (60, '("            - hmllolon", T40,": ", F13.2)') time_hmllolon
        Write (60, '("            - olpaan", T40,": ", F13.2)') time_olpaan
        Write (60, '("            - olpalon", T40,": ", F13.2)') time_olpalon
        Write (60, '("            - olplolon", T40,": ", F13.2)') time_olplolon
        Write (60, '("            - hmlistln", T40,": ", F13.2)') time_hmlistln
        Write (60, '("            - olpistln", T40,": ", F13.2)') time_olpistln

        Write (60, '(" first-variational secular equation", T40, ": ", F12.2)') timefv
        If (associated(input%groundstate%spin)) Then
          Write (60, '(" second-variational calculation", T40, ": ", F12.2)') timesv
        End If
        Write (60, '(" charge density calculation", T40, ": ", F12.2)') timerho
        Write (60, '(" potential calculation", T40, ": ", F12.2)') timepot
        Write (60, '(" muffin-tin manipulations", T40, ": ", F12.2)') timemt
        Write (60, '(" APW matching", T40, ": ", F12.2)') timematch
        Write (60, '(" disk reads/writes", T40, ": ", F12.2)') timeio
        Write (60, '(" mixing efforts", T40, ": ", F12.2)') timemixer
        If (input%groundstate%tforce) Then
          Write (60, '(" force calculation", T40, ": ", F12.2)') timefor
        End If
        timetot = timeinit + timemat + timefv + timesv + timerho + &
       &          timepot + timefor + timeio + timemt + timemixer + timematch
        Write (60, '(" sum", T40, ": ", F12.2)') timetot
        Call timesec(tsg1)
        Write (60, '(" total", T40, ": ", F12.2)') tsg1-tsg0
        Write (60,*)
        Write (60, '("More timings (seconds)!")')
        Write (60, '(" Dirac eqn solver", T40, ": ", F12.2)') time_rdirac
        Write (60, '(" Rel. Schroedinger eqn solver", T40, ": ", F12.2)') time_rschrod
        Write (60, '(" Total time spent in radial solvers", T40, ": ", F12.2)') time_rdirac+time_rschrod
        Write (60, '("+----------------------------+")')
        Write (60, '("| Groundstate module stopped |")')
        Write (60, '("+----------------------------+")')
        close (60)
    endif

    Return
End Subroutine gndstate
!EOC
