!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine phdos
      Use modmain
      Use modinput
      Use FoX_wxml
      Use modmpi
      Use physical_constants, only: kboltz
      use phonons_util, only: ph_util_setup_interpolation, ph_util_interpolate, ph_util_diag_dynmat
      use phonons_io_util, only: ph_io_read_dielten, ph_io_read_borncharge
      use matrix_fourier_interpolation, only: mfi_type
      use mod_kpointset, only: k_set, generate_k_vectors, delete_k_vectors
      use mod_opt_tetra
      Implicit None
! local variables
      Integer :: n, iq, i, iw, ngridqint(3)
      Real (8) :: wmin, wmax, wd, dw, dielten(3, 3)
      Real (8) :: tmax, temp (input%phonons%phonondos%ntemp), s (input%phonons%phonondos%ntemp)
      Real (8) :: t1, t2
      type(k_set) :: qset_int
      type(t_set) :: tset
      type(mfi_type) :: mfi
      Type (xmlf_t), Save :: xf
      Character (256) :: buffer
      logical :: success
! allocatable arrays
      Real (8), Allocatable :: borncharge(:, :, :)
      Real (8), Allocatable :: wp (:, :)
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: wgt (:, :, :)
      Real (8), Allocatable :: gw (:)
      Real (8), Allocatable :: f (:), g (:), cf (:, :)
      Complex (8), Allocatable :: dynr (:, :, :)
      Complex (8), Allocatable :: dynp (:, :, :)
      Complex (8), Allocatable :: ev (:, :)
      ! do this only in master process
      if (rank .ne. 0) goto 10
      ! initialise universal variables
      Call init0
      Call init2
      n = 3 * natmtot
      ! try to read dielectric tensor and Born effective charges
      if( input%phonons%polar ) then
        call ph_io_read_dielten( dielten, 'EPSINF.OUT', success )
        call terminate_if_false( success, '(phdos) &
          Failed to read dielectric tensor from EPSINF.OUT.' )
        call ph_io_read_borncharge( borncharge, 'ZSTAR.OUT', success )
        call terminate_if_false( success, '(phdos) &
          Failed to read Born effective charge tensors from ZSTAR.OUT.' )
      end if
      ! set up interpolation
      if( input%phonons%polar ) then
        call ph_util_setup_interpolation( bvec, ngridq, nqpt, ivq(:, 1:nqpt), vql(:, 1:nqpt), mfi, dynr, &
          dielten=dielten, borncharge=borncharge, &
          sumrule=input%phonons%sumrule )
      else
        call ph_util_setup_interpolation( bvec, ngridq, nqpt, ivq(:, 1:nqpt), vql(:, 1:nqpt), mfi, dynr, &
          sumrule=input%phonons%sumrule )
      end if
      ! generate integration q-grid
      ngridqint = input%phonons%phonondos%ngridqint
      if( input%phonons%phonondos%ngrdos > 0 ) &
        ngridqint = input%phonons%phonondos%ngrdos
      call generate_k_vectors( qset_int, bvec, ngridqint, [0.d0, 0.d0, 0.d0], .true., uselibzint=.false. )
      ! interpolate phonon frequencies on integration grid
      if( input%phonons%polar ) then
        call ph_util_interpolate( qset_int%nkpt, qset_int%vkl(:, 1:qset_int%nkpt), mfi, dynr, dynp, &
          dielten=dielten, borncharge=borncharge )
      else
        call ph_util_interpolate( qset_int%nkpt, qset_int%vkl(:, 1:qset_int%nkpt), mfi, dynr, dynp )
      end if
      allocate( wp(n, qset_int%nkpt), ev(n, n) )
      do iq = 1, qset_int%nkpt
        call ph_util_diag_dynmat( dynp(:, :, iq), wp(:, iq), ev )
      end do
      ! find the minimum and maximum frequencies
      wmin = minval( wp )
      wmax = maxval( wp )
      wmax = wmax + (wmax-wmin) * 0.1d0
      wmin = wmin - (wmax-wmin) * 0.1d0
      ! generate energy grid
      wd = wmax - wmin
      if( wd < 1.d-8 ) wd = 1.d0
      dw = wd / input%phonons%phonondos%nwdos
      allocate( w(input%phonons%phonondos%nwdos) )
      do iw = 1, input%phonons%phonondos%nwdos
        w(iw) = dw * dble(iw-1) + wmin
      end do
      ! calculate DOS
      allocate( gw(input%phonons%phonondos%nwdos), source=0.d0 )
      ! new method (improved tetrahedron integration)
      if( input%phonons%phonondos%inttype == 'tetra' ) then
        call opt_tetra_init( tset, qset_int, 2, reduce=.true. )
        allocate( wgt( n, qset_int%nkpt, input%phonons%phonondos%nwdos) )
        call opt_tetra_wgt_delta( tset, qset_int%nkpt, n, wp, input%phonons%phonondos%nwdos, w, wgt )
        do iw = 1, input%phonons%phonondos%nwdos
          gw(iw) = sum( wgt(:, : ,iw) )
        end do
        deallocate( wgt )
        call opt_tetra_destroy( tset )
      ! old method (direct summation)
      else
        !$omp parallel default( shared ) private( t1, iw ) reduction( +: gw )
        !$omp do collapse(2)
        do iq = 1, qset_int%nkpt
          do i = 1, n
            t1 = (wp(i, iq) - wmin) / dw + 1.d0
            iw = nint( t1 )
            if( (iw >= 1) .and. (iw <= input%phonons%phonondos%nwdos) ) &
              gw(iw) = gw(iw) + qset_int%wkpt(iq)
          end do
        end do
        !$omp end do
        !$omp end parallel
        t1 = 1.d0 / dw
        gw = t1 * gw
      end if
      ! smooth phonon DOS if required
      if( input%phonons%phonondos%nsmdos > 0) &
        call fsmooth( input%phonons%phonondos%nsmdos, input%phonons%phonondos%nwdos, 1, gw )
      ! write phonon DOS to file
      open( 50, file='PHDOS.OUT', action='WRITE', form='FORMATTED' )
      do iw = 1, input%phonons%phonondos%nwdos
        write( 50, '(2G18.10)' ) w(iw), gw(iw)
      end do
      close( 50 )
      write( *, * )
      write( *, '("Info(phdos):")' )
      write( *, '(" phonon density of states written to PHDOS.OUT")' )
      ! clean up
      deallocate( wp, ev )
      if( allocated( dynr ) ) deallocate( dynr )
      if( allocated( dynp ) ) deallocate( dynp )
      if( allocated( borncharge ) ) deallocate( borncharge )
      call delete_k_vectors( qset_int )
      call mfi%destroy
!-------------------------------------------!
!     thermodynamic properties from DOS     !
!-------------------------------------------!
      Allocate (f(input%phonons%phonondos%nwdos), g(input%phonons%phonondos%nwdos), cf(3, input%phonons%phonondos%nwdos))
! maximum temperature
      tmax = wmax / kboltz
! temperature grid
      Do i = 1, input%phonons%phonondos%ntemp
         temp (i) = tmax * dble (i) / dble (input%phonons%phonondos%ntemp)
      End Do
      Open (50, File='THERMO.OUT', Action='WRITE', Form='FORMATTED')
      Call xml_OpenFile ("thermo.xml", xf, replace=.True., pretty_print=.True.)
      Call xml_NewElement (xf, "thermodynamicproperties")
! zero point energy
      Do iw = 1, input%phonons%phonondos%nwdos
         f (iw) = gw (iw) * w (iw)
      End Do
      Call fderiv (-1, input%phonons%phonondos%nwdos, w, f, g, cf)
      t1 = 0.5d0 * g (input%phonons%phonondos%nwdos)
      Write (50,*)
      Write (50, '("Zero-point energy : ", G18.10)') t1
      Call xml_NewElement (xf, "zeropointenergy")
      Write (buffer, '(g18.10)') t1
      Call xml_AddAttribute (xf, "name", "zero-point energy")
      Call xml_AddAttribute (xf, "value", trim(adjustl(buffer)))
      Call xml_AddAttribute (xf, "unit", "Hartree")
      Call xml_endElement (xf, "zeropointenergy")
! vibrational energy
      Write (50,*)
      Write (50, '("Vibrational energy vs. temperature :")')
      Call xml_NewElement (xf, "vibrationalenergy")
      Call xml_NewElement (xf, "mapdef")
      Call xml_NewElement (xf, "variable1")
      Call xml_AddAttribute (xf, "name", "temperature")
      Call xml_AddAttribute (xf, "unit", "Kelvin")
      Call xml_endElement (xf, "variable1")
      Call xml_NewElement (xf, "function1")
      Call xml_AddAttribute (xf, "name", "vibrational energy")
      Call xml_AddAttribute (xf, "unit", "Hartree")
      Call xml_endElement (xf, "function1")
      Call xml_endElement (xf, "mapdef")
      Do i = 1, input%phonons%phonondos%ntemp
         Do iw = 1, input%phonons%phonondos%nwdos
            t1 = w (iw) / (2.d0*kboltz*temp(i))
            If (t1 .Gt. 0.d0) Then
               f (iw) = gw (iw) * w (iw) * Cosh (t1) / Sinh (t1)
            Else
               f (iw) = 0.d0
            End If
         End Do
         Call fderiv (-1, input%phonons%phonondos%nwdos, w, f, g, cf)
         t1 = 0.5d0 * g (input%phonons%phonondos%nwdos)
         Write (50, '(2G18.10)') temp (i), t1
         Call xml_NewElement (xf, "map")
         Write (buffer, '(4g18.10)') temp(i)
         Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
         Write (buffer, '(4g18.10)') t1
         Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
         Call xml_endElement (xf, "map")
         s (i) = t1
      End Do
      Call xml_endElement (xf, "vibrationalenergy")
! free energy
      Write (50,*)
      Write (50, '("Free energy vs. temperature :")')
      Call xml_NewElement (xf, "vibrationalfreeenergy")
      Call xml_NewElement (xf, "mapdef")
      Call xml_NewElement (xf, "variable1")
      Call xml_AddAttribute (xf, "name", "temperature")
      Call xml_AddAttribute (xf, "unit", "Kelvin")
      Call xml_endElement (xf, "variable1")
      Call xml_NewElement (xf, "function1")
      Call xml_AddAttribute (xf, "name", "vibrational free energy")
      Call xml_AddAttribute (xf, "unit", "Hartree")
      Call xml_endElement (xf, "function1")
      Call xml_endElement (xf, "mapdef")
      Do i = 1, input%phonons%phonondos%ntemp
         Do iw = 1, input%phonons%phonondos%nwdos
            t1 = 2.d0 * Sinh (w(iw)/(2.d0*kboltz*temp(i)))
            If (t1 .Gt. 0.d0) Then
               f (iw) = gw (iw) * Log (t1)
            Else
               f (iw) = 0.d0
            End If
         End Do
         Call fderiv (-1, input%phonons%phonondos%nwdos, w, f, g, cf)
         t1 = kboltz * temp (i) * g &
        & (input%phonons%phonondos%nwdos)
         Write (50, '(2G18.10)') temp (i), t1
         Call xml_NewElement (xf, "map")
         Write (buffer, '(4g18.10)') temp(i)
         Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
         Write (buffer, '(4g18.10)') t1
         Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
         Call xml_endElement (xf, "map")
! compute entropy from S = (F-E)/T
         s (i) = (s(i)-t1) / temp (i)
      End Do
      Call xml_endElement (xf, "vibrationalfreeenergy")
! entropy
      Write (50,*)
      Write (50, '("Entropy vs. temperature :")')
      Call xml_NewElement (xf, "vibrationalentropy")
      Call xml_NewElement (xf, "mapdef")
      Call xml_NewElement (xf, "variable1")
      Call xml_AddAttribute (xf, "name", "temperature")
      Call xml_AddAttribute (xf, "unit", "Kelvin")
      Call xml_endElement (xf, "variable1")
      Call xml_NewElement (xf, "function1")
      Call xml_AddAttribute (xf, "name", "vibrational entropy")
      Call xml_AddAttribute (xf, "unit", "Hartree/Kelvin")
      Call xml_endElement (xf, "function1")
      Call xml_endElement (xf, "mapdef")
      Do i = 1, input%phonons%phonondos%ntemp
         Write (50, '(2G18.10)') temp (i), s (i)
         Call xml_NewElement (xf, "map")
         Write (buffer, '(4g18.10)') temp(i)
         Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
         Write (buffer, '(4g18.10)') s(i)
         Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
         Call xml_endElement (xf, "map")
      End Do
      Call xml_endElement (xf, "vibrationalentropy")
! heat capacity
      Write (50,*)
      Write (50, '("Heat capacity vs. temperature :")')
      Call xml_NewElement (xf, "heatcapacity")
      Call xml_NewElement (xf, "mapdef")
      Call xml_NewElement (xf, "variable1")
      Call xml_AddAttribute (xf, "name", "temperature")
      Call xml_AddAttribute (xf, "unit", "Kelvin")
      Call xml_endElement (xf, "variable1")
      Call xml_NewElement (xf, "function1")
      Call xml_AddAttribute (xf, "name", "heat capacity")
      Call xml_AddAttribute (xf, "unit", "Hartree/Kelvin")
      Call xml_endElement (xf, "function1")
      Call xml_endElement (xf, "mapdef")
      Do i = 1, input%phonons%phonondos%ntemp
         Do iw = 1, input%phonons%phonondos%nwdos
            t1 = w (iw) / (kboltz*temp(i))
            t2 = Exp (t1) - 1.d0
            If (t2 .Ne. 0.d0) Then
               f (iw) = gw (iw) * (t1**2) * (t2+1.d0) / t2 ** 2
            Else
               f (iw) = 0.d0
            End If
         End Do
         Call fderiv (-1, input%phonons%phonondos%nwdos, w, f, g, cf)
         t1 = kboltz * g (input%phonons%phonondos%nwdos)
         Write (50, '(2G18.10)') temp (i), t1
         Call xml_NewElement (xf, "map")
         Write (buffer, '(4g18.10)') temp(i)
         Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
         Write (buffer, '(4g18.10)') t1
         Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
         Call xml_endElement (xf, "map")
      End Do
      Call xml_endElement (xf, "heatcapacity")
      Close (50)
      Call xml_endElement (xf, "thermodynamicproperties")
      Call xml_Close (xf)
      Write (*, '(" thermodynamic properties written to THERMO.OUT")')
      Write (*,*)
      Deallocate ( w, gw, f, g, cf )
10    continue
#ifdef MPI
      call MPI_Barrier(MPI_Comm_World, ierr)
#endif
      Return
End Subroutine
