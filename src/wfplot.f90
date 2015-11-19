!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine wfplot(dostm)
    Use modmain
    Use modinput
    use modplotlabels
    use modmpi, only : rank
    Implicit None
    Logical, Intent(in) :: dostm
    ! local variables
    Integer :: ik, ist, stmmode, stmtype
    Real (8) :: x, y, t1, bias
    type(plotlabels),pointer ::labels
    ! allocatable arrays
    Complex (8), Allocatable :: evecfv (:, :)
    Complex (8), Allocatable :: evecsv (:, :)
    character(256) :: string
    ! external functions
    Real (8) :: sdelta, stheta
    External sdelta, stheta
    ! initialise universal variables
    Call init0
    Call init1
    Allocate (evecfv(nmatmax, nstfv))
    Allocate (evecsv(nstsv, nstsv))
! initialise the charge density and potentials from file
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
! in case of HF hybrids use PBE potential
            string=filext
            filext='_PBE.OUT'
            Call readstate
            filext=string
           Else
               Call readstate
           End If
        Else         
           Call readstate
        End If 
! read Fermi energy from file
    Call readfermi
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! update potential in case if HF Hybrids
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
               Call readstate
           End If
        End If 
    ! set the occupancies
    If ( .Not. dostm) Then
        ! kstlist should only contain one k-point and state for wave-function plot
        if (size(input%properties%wfplot%kstlist%pointstatepair,2).ne.1) then
          if (rank==0) then
            write(*,*)
            write(*,'("Error(wfplot): /input/properties/wfplot/kstlist must contain")')
            write(*,'(" only one pointstatepair, but ",i6," were defined")') &
                size(input%properties%wfplot%kstlist%pointstatepair,2)
            write(*,*)
            stop
          end if
        end if
        ik = input%properties%wfplot%kstlist%pointstatepair(1, 1)
        ist = input%properties%wfplot%kstlist%pointstatepair(2, 1)
        If ((ik .Lt. 1) .Or. (ik .Gt. nkpt)) Then
          if (rank==0) then
            Write (*,*)
            Write (*, '("Error(wfplot): k-point out of range : ", I8)') &
                & ik
            Write (*,*)
            Stop
          end if
        End If
        If ((ist .Lt. 1) .Or. (ist .Gt. nstsv)) Then
          if (rank==0) then
            Write (*,*)
            Write (*, '("Error(wfplot): state out of range : ", I8)') &
                & ist
            Write (*,*)
            Stop
          end if
        End If
        ! plotting a single wavefunction
        occsv (:, :) = 0.d0
        occsv (ist, ik) = 1.d0
    Else
        Select Case (trim(input%properties%STM%stmtype))
            Case ('differentialConductance')
                stmtype=1
            Case ('integratedLDOS')
                stmtype=2
        End Select
        Select Case (trim(input%properties%STM%stmmode))
            Case ('constantHeight')
                stmmode=1
            Case ('topographic')
                stmmode=2
        End Select

        bias = input%properties%STM%bias

        If ( stmtype .Eq. 1 .And. stmmode .Eq. 1) Then
          if (rank==0) then
            Write(*,*)
            Write (*, '("Info(wfplot):")')
            Write (*, '("Generating constant-height STM image of the differential conductance.")')
            ! plotting an STM differential-conductance image by setting occupancies to be a
            ! delta function at the Fermi energy
            t1 = 1.d0 / input%groundstate%swidth
            Do ik = 1, nkpt
                ! get the eigenvalues from file
                Call getevalsv (vkl(:, ik), evalsv(:, ik))
                Do ist = 1, nstsv
                    x = ((efermi-bias)-evalsv(ist, ik)) * t1
                    occsv (ist, ik) = occmax * wkpt (ik) * sdelta &
                        & (input%groundstate%stypenumber, x) * t1
                End Do
            End Do
          end if
        Else If ( stmtype .Eq. 2 .And. stmmode .Eq. 1) Then
           ! Plots the local density of states integrated between Ef y Ef + bias for positive bias or
           ! between Ef-bias and Ef for negative bias. This way simple STM plot in the Tersoff-Hamann
           ! aproximation can be obtained (PRB 31,805 (1985)).
          if (rank==0) then
            Write(*,*)
            Write (*, '("Info(wfplot):")')
            Write (*, '("Generating constant-height STM image of the integrated LDOS.")')
            t1 = 1.d0 / input%groundstate%swidth
            Do ik = 1, nkpt
                Call getevalsv (vkl(:, ik), evalsv(:, ik))
                Do ist = 1, nstsv
                    x = sign(1.d0,bias)*(efermi+bias-evalsv(ist, ik)) * t1
                    y = sign(1.d0,bias)*(evalsv(ist, ik)-efermi) * t1
                    occsv (ist, ik) = occmax * wkpt (ik) * stheta &
                        & (input%groundstate%stypenumber, x)*&
                        & stheta(input%groundstate%stypenumber, y)
                End Do
            End Do
          end if
        Else
            call warning('Warning(wfplot): STM still not implemented for direct topographic plot.')
            call warning('For topographic plot generation consider to make a series &
              constant-height calculations at different heights and postprocess the &
              output to find the iso-surface.')
        End If
    End If
    ! set the charge density to zero
    rhomt (:, :, :) = 0.d0
    rhoir (:) = 0.d0
    ! compute the charge density with the new occupancies
    Do ik = 1, nkpt
        ! get the eigenvectors from file
        Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
        Call getevecsv (vkl(:, ik), evecsv)
        Call rhovalk (ik, evecfv, evecsv)
        Call genrhoir (ik, evecfv, evecsv)
    End Do
    ! symmetrise the density for the STM plot
    If (dostm) Then
        Call symrf (input%groundstate%lradstep, rhomt, rhoir)
    End If
    ! convert the density from a coarse to a fine radial mesh
    Call rfmtctof (rhomt)

    ! write the wavefunction modulus squared plot to file
    If (associated(input%properties%wfplot)) Then

        If (associated(input%properties%wfplot%plot1d)) Then
            labels=>create_plotlablels("Potential","WF1D",1)
            call set_plotlabel_axis(labels,1,"Distance","a_0","graceunit")
            call set_plotlabel_axis(labels,2,"Probability Density","","graceunit")
            Call plot1d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
                & rhomt, rhoir, input%properties%wfplot%plot1d)
            call destroy_plotlablels(labels)
            if (rank==0) then
              Write (*,*)
              Write (*, '("Info(wfplot):")')
              Write (*, '(" 1D wavefunction modulus squared written to WF1D.xml")')
            end if
        End If

        If (associated(input%properties%wfplot%plot2d)) Then
            labels=>create_plotlablels("Probability Density","WF2D",2)
            call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
            call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
            call set_plotlabel_axis(labels,3,"Wave Function Norm Squared","","graceunit")
            Call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
                & rhomt, rhoir, input%properties%wfplot%plot2d)
            call destroy_plotlablels(labels)
            if (rank==0) then
              Write (*,*)
              Write (*, '("Info(wfplot):")')
              Write (*, '(" 2D wavefunction modulus squared written to WF2D.xml")')
            end if
        End If

        If (associated(input%properties%wfplot%plot3d)) Then
            input%properties%wfplot%plot3d%usesym=.false.
            labels=>create_plotlablels("Probability Density","WF3D",3)
            call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
            call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
            call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
            call set_plotlabel_axis(labels,4,"Wave Function Norm Squared","","graceunit")
            Call plot3d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
                & rhomt, rhoir, input%properties%wfplot%plot3d)
            call destroy_plotlablels(labels)
            if (rank==0) then
              Write(*,*)
              Write(*, '("Info(wfplot):")')
              Write(*, '(" 3D wavefunction modulus squared written to WF3D.xml")')
            end if
        End If

    End If

    If (dostm) Then
        labels=>create_plotlablels("2D STM image","STM2d",2)
        call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
        call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
        call set_plotlabel_axis(labels,3,"STM","","graceunit")
        Call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
            & rhomt, rhoir, input%properties%STM%plot2d)
        call destroy_plotlablels(labels)
        If ( stmmode .Eq. 2) Then
          if (rank==0) then
            Write (*,*)
            Write (*, '("Info(wfplot):")')
            Write (*, '("STM still not implemented for direct topographic plot.")')
            Write (*, '("For topographic plot generation consider to make a series &
            constant-height calculations at different heights and postprocess the &
               output to find the iso-surface. ")')
          end if
        Else
          if (rank==0) then
            Write (*,*)
            Write (*, '("Info(wfplot):")')
            Write (*, '(" 2D STM image written to STM2d2D.xml")')
          end if
        End If
    End If
    If ( .Not. dostm) Then
      if (rank==0) then
        Write (*, '(" for k-point ", I6, " and state ", I6)') &
            input%properties%wfplot%kstlist%pointstatepair(1, 1), &
            input%properties%wfplot%kstlist%pointstatepair(2, 1)
      end if
    End If
    if (rank==0) Write (*,*)
    Deallocate (evecfv, evecsv)
    Return
End Subroutine

