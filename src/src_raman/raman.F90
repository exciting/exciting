! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted for exciting 2014, Stefan Kontur
!
! =============================================================================
subroutine RAMAN
! =============================================================================
!
!     determines eigenvalues for an arbitrary potential
!     V(x)= A0 + A1 x + A2 x**2 + A3 x**3 + A4 x**4 + A5 x**5 + A6 x**6 
!     and computes Raman scattering intensities
!
use mod_lattice, only: omega
use mod_misc
use mod_qpoint, only: nqpt, ngridq
use mod_atoms, only: natmtot, nspecies, natoms, spmass, idxas
use mod_energy, only: engytot
use mod_force, only: forcetot
use modinput
use modxs
use raman_coeff
use raman_ew
use raman_inter
use raman_fij
use raman_trmat
use m_symvec
use m_raman_utils
use m_genfilname
#ifdef IFORT
  use ifport
#endif
#ifdef MPI
  use modmpi
#endif
!
implicit none
! Raman data
integer :: maxp
integer :: i, j, ia, is, iat, ic, imode, istep, nmode, iw
integer :: istep_lo, istep_hi, i_shift
integer :: oct1, oct2, ias
integer :: read_i
real(8) :: rlas,sn,zmin,zs, norm
real(8) :: dph, force_sum, vgamc(3)
real(8) :: read_dph, read_engy
real(8) :: wlas, ws, dwlas, Sab, t1
! timing
real(8) :: start_time_cpu, finish_time_cpu, time_cpu_tot, t_cpu_proc
real(8) :: start_time_wall, finish_time_wall, time_wall_tot
! phonons
Real(8), Allocatable :: w(:)
Complex(8), Allocatable :: ev(:, :)
Complex(8), Allocatable :: dynq(:, :, :)
Complex(8), Allocatable :: dynp(:, :)
Complex(8), Allocatable :: dynr(:, :, :)
! check modes
Logical, Allocatable :: active(:), acoustic(:)
integer, allocatable :: irep(:)
Logical :: existent, existent1, eq_done, nlf
! file handling
Character(256) :: raman_filext, raman_stepdir, fileeps(3, 3)
! check whether plan was specified in xs element
Logical :: dplan
Integer :: nxstasksmax, taskindex
! BSE variables
Logical :: fcoup
character(256) :: frmt, tdastring, bsetypestring, scrtypestring, &
  &               epsilondir
#ifndef IFORT
  integer system, chdir
#endif
! MPI
#ifndef MPI
  integer, parameter :: rank = 0
#endif
!
!
! we require an XS input
if (.not. associated(input%xs)) then
   write(*, '("Error(Raman): Raman runs require the specification of the input%xs element!")')
   stop
endif
!
! take time
time_cpu_tot = 0.d0; time_wall_tot = 0.d0
call cpu_time(t_cpu_proc)
#ifdef MPI
   call MPI_Allreduce(t_cpu_proc, start_time_cpu, 1, MPI_Real8, MPI_Sum, MPI_Comm_World, ierr)
#else
   start_time_cpu = t_cpu_proc
#endif
call timesec(start_time_wall)
!
call init0
!
! default number of modes
nmode = 3*natmtot
!
!
! open INFO file for all pre-Raman computations triggered in this subroutine
if (rank .eq. 0) then
   open(unit=99,file='INFO_RAMAN.OUT',status='unknown',form='formatted')
endif
!
!
numpot = input%properties%raman%nstep
allocate( potinx(3*natmtot, numpot) )
allocate( potiny(3*natmtot, numpot) )
!
! scattering volume for molecules or solids
if (input%properties%raman%molecule) then
   ncell = 1
else
   ncell = 100000
endif
!
!     laser energy   
select case(trim(input%properties%raman%elaserunit))
  case ('eV')
     rlas = input%properties%raman%elaser * fevha
  case ('nm')
     if (input%properties%raman%elaser .gt. 1.0d-6) then
        rlas = 1.d0/input%properties%raman%elaser * frnmha
     else 
        rlas = 1.d6 * frnmha
     endif
  case ('cm-1')
     rlas = input%properties%raman%elaser * fwnha
  case default
end select
!     rlas is now in Ha
!
!
!
nwdf = input%xs%energywindow%points
if (associated (input%xs%tddft)) then
   If (input%xs%tddft%acont) Then
      nwdf = input%xs%tddft%nwacont
   End If
endif
!
!
!
! == allocate arrays
!
! eigenvalues and vectors, dynamical matrices
Allocate (w(3*natmtot))
Allocate (ev(3*natmtot, 3*natmtot))
Allocate (active(3*natmtot))
Allocate (acoustic(3*natmtot))
Allocate (irep(3*natmtot))
! dielectric functions
allocate( df(3*natmtot, numpot, 3, 3, nwdf) )
! arrays internal to Raman computation
allocate( eigen(2*input%properties%raman%ninter) )
allocate( T(2*input%properties%raman%ninter,2*input%properties%raman%ninter) )
allocate( R(2*input%properties%raman%ninter,2*input%properties%raman%ninter) )
allocate( z1(input%properties%raman%ninter+1,input%properties%raman%nstate) )
allocate( z2(input%properties%raman%ninter+1,input%properties%raman%nstate) )
allocate( transme1(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( transme2(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( transme3(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( transme4(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( transme5(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( transme6(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( e1(input%properties%raman%nstate) )
allocate( e2(input%properties%raman%nstate) )
allocate( e3(input%properties%raman%nstate) )
allocate( de(input%properties%raman%nstate) )
allocate( indexi(input%properties%raman%nstate) )
allocate( xa(input%properties%raman%ninter) )
allocate( xpot(4*input%properties%raman%ninter+1) )
allocate( pot(4*input%properties%raman%ninter+1) )
allocate( b0(input%properties%raman%ninter) )
allocate( b1(input%properties%raman%ninter) )
allocate( b2(input%properties%raman%ninter) )
allocate( b3(input%properties%raman%ninter) )
allocate( b4(input%properties%raman%ninter) )
!
!check whether plan element is specified in xs element
dplan=associated(input%xs%plan)
if (rank .eq. 0) then
   ! header for INFO_RAMAN.OUT
   Write (99, '("+-----------------------------------------------------------+")')
   Write (99, '("| EXCITING ",A)') versionname
   Write (99, '("| version hash id: ",a)') githash
#ifdef MPI
   Write (99, '("| MPI version using ",i6," processor(s)                     |")') procs
#ifndef MPI1
   Write (99, '("|  using MPI-2 features                                     |")')
#endif
#endif
   Write (99, '("+-----------------------------------------------------------+",//)')
   Write (99, '("+-----------------------------------------------------------+")')
   Write (99, '("|       Info file for all computation steps necessary       |")')
   Write (99, '("|             in the workflow for calculating               |")')
   Write (99, '("|              Raman scattering intensities                 |")')
   Write (99, '("+-----------------------------------------------------------+",//)')
endif
!
! save input structure
call raman_save_struct
!
! +++ SYMMETRY ANALYSIS OF THE CRYSTAL +++
! +++ CONSTRUCTION OF CHARACTER TABLE  +++
! 
input%structure%tshift = .true.
if (rank .eq. 0) then
   call construct_chartabl
   write(99,'("Info(Raman): Character table constructed.",/)')
endif
!
! +++ TRIGGER ALL NECESSARY COMPUTATIONS FROM OTHER EXCITING PARTS TO +++
! +++     COLLECT DATA FOR POTENTIALS AND DIELECTRIC FUNCTIONS        +++
!
! get normal coordinates for phonons at Gamma
!
vgamc(:) = 0.d0
!
Select Case (input%properties%raman%getphonon)
   Case ('fromscratch')
      ! do supercell phonon calculation
      if (rank .eq. 0) then
         write(99,'("Info(Raman): Performing a supercell phonon calculation for the Gamma point.",/)')
         call flushifc(99)
      endif
      ! create dummy phonon element if not present in the input file
      If ( .Not. (associated(input%phonons))) Then
         input%phonons => getstructphonons (emptynode)
      End If
      If (associated(input%phonons%parts)) nullify (input%phonons%parts)
      ! set parameters for Gamma phonons
      input%phonons%do = 'fromscratch'
      input%phonons%ngridq = (/ 1, 1, 1/)
      task = 200
      call phonon
      ! restore original structure
      call raman_restore_struct
      ! continue with reading from files... (q-points were initialized in phonon)
      Allocate (dynq(3*natmtot, 3*natmtot, nqpt))
      Allocate (dynp(3*natmtot, 3*natmtot))
      Allocate (dynr(3*natmtot, 3*natmtot, ngridq(1)*ngridq(2)*ngridq(3)))
      Call readdyn (.true.,dynq)
      ! apply the acoustic sum rule
      Call sumrule (dynq)
      ! obtain eigenvectors at Gamma
      Call dynqtor (dynq, dynr)
      Call dynrtoq (vgamc(:), dynr, dynp)
      Call dyndiag (dynp, w, ev)
      if (rank .eq. 0) then
       do imode = 1, 3*natmtot
          Write(99, '(/," Mode ",i3)') imode
          write(99, '(" Eigenvalue: ",f15.5," cm^-1",/)') w(imode)*fhawn
          Write(99, '("      Atom   Polarization    Eigenvector")')
          do iat = 1, natmtot
             do i = 1, 3
                write(99, '(i10,i15,2f12.7)') iat, i, dble(ev(3*(iat-1)+i, imode)), aimag(ev(3*(iat-1)+i, imode))
             enddo
          enddo
       enddo
       Write(99, *)
       call flushifc(99)
      endif
      ! reset file extension to default
      filext = '.OUT'
!
   Case ('fromfile')
   ! read in the dynamical matrix from files
   if (rank .eq. 0) write(99,'("Info(Raman): Reading dynamical matrix from files.",/)')
      ! create dummy phonon element if not present in the input file
      If ( .Not. (associated(input%phonons))) Then
         input%phonons => getstructphonons (emptynode)
      End If
      ! initialize q-dependent variables
      task = 200
      call init2
      Allocate (dynq(3*natmtot, 3*natmtot, nqpt))
      Allocate (dynp(3*natmtot, 3*natmtot))
      Allocate (dynr(3*natmtot, 3*natmtot, ngridq(1)*ngridq(2)*ngridq(3)))
      ! read dynamical matrix
      Call readdyn (.true.,dynq)
      ! apply the acoustic sum rule
      Call sumrule (dynq)
      ! obtain eigenvectors at Gamma
      Call dynqtor (dynq, dynr)
      Call dynrtoq (vgamc(:), dynr, dynp)
      Call dyndiag (dynp, w, ev)
      if (rank .eq. 0) then
       do imode = 1, 3*natmtot
          Write(99, '(/," Mode ",i3)') imode
          write(99, '(" Eigenvalue: ",f15.5," cm^-1",/)') w(imode)*fhawn
          Write(99, '("      Atom   Polarization    Eigenvector")')
          do iat = 1, natmtot
             do i = 1, 3
                write(99, '(i10,i15,2f12.7)') iat, i, dble(ev(3*(iat-1)+i, imode)), aimag(ev(3*(iat-1)+i, imode))
             enddo
          enddo
       enddo
       Write(99, *)
       call flushifc(99)
      endif
!
   Case ('symvec')
   ! construct symmetry vectors and use them as normal coordinates
      if (rank .eq. 0) then
       call construct_symvec(ev)
       write(99,'("Info(Raman): Symmetry vectors constructed:",/)')
       do imode = 1, 3*natmtot
          Write(99, '(/," Mode ",i3)') imode
          Write(99, '("      Atom   Polarization    Eigenvector")')
          do iat = 1, natmtot
             do i = 1, 3
                write(99, '(i10,i15,2f12.7)') iat, i, dble(ev(3*(iat-1)+i, imode)), aimag(ev(3*(iat-1)+i, imode))
             enddo
          enddo
       enddo
       Write(99, *)
       call flushifc(99)
      endif
#ifdef MPI
      ! broadcast eigenvectors to child processes
      call MPI_Bcast(ev, 3*natmtot*3*natmtot, MPI_Double_Complex, 0, MPI_Comm_World, ierr)
#endif
!
   Case ('symveccheck')
   ! construct symmetry vectors and return
      if (rank .eq. 0) then
       call construct_symvec(ev)
       write(99,'("Info(Raman): Symmetry vectors constructed:",/)')
       do imode = 1, 3*natmtot
          Write(99, '(/," Mode ",i3)') imode
          Write(99, '("      Atom   Polarization    Eigenvector")')
          do iat = 1, natmtot
             do i = 1, 3
                write(99, '(i10,i15,2f12.7)') iat, i, dble(ev(3*(iat-1)+i, imode)), aimag(ev(3*(iat-1)+i, imode))
             enddo
          enddo
       enddo
       Write(99, '(/,"Info(Raman): getphonon=symveccheck, so we stop here ",/)')
      endif
      if (associated(input%phonons)) nullify (input%phonons)
      if (associated(input%gw)) nullify (input%gw)
      nullify (input%xs)
      return
   Case ('readinput')
      if (rank .eq. 0) write(99,'("Info(Raman): Phonon eigenvectors read from input file.",/)')
      ! read from input file
      ic = 0
      norm = 0.d0
      ev = zzero
      do iat = 1, natmtot
         do i = 1, 3
            ic = ic + 1
            ev(ic, 1) = cmplx(input%properties%raman%eigvecarray(ic)%eigvec%comp(1), &
                  &           input%properties%raman%eigvecarray(ic)%eigvec%comp(2), 8)
            norm = norm + abs(ev(ic, 1))**2
         enddo
      enddo
      if (abs(norm - 1.d0) .gt. eps .and. norm .gt. eps) then
       if (rank .eq. 0) then
         write(99, '("Info(Raman): Norm of eigenvector is ",f14.9)') sqrt(norm)
         write(99, '("             Eigenvector will be renormalized")')
         call flushifc(99)
       endif
       ev = ev / sqrt(norm)
      endif
      ! use input variable mode for direct loop over modes
      nmode = 1
      input%properties%raman%mode = 0
End Select
!
! take time
call cpu_time(t_cpu_proc)
#ifdef MPI
  call MPI_Allreduce(t_cpu_proc, finish_time_cpu, 1, MPI_Real8, MPI_Sum, MPI_Comm_World, ierr)
#else 
  finish_time_cpu = t_cpu_proc
#endif
call timesec(finish_time_wall)
time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
if (rank .eq. 0) then
  write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') time_cpu_tot, time_cpu_tot/3600.
  write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') time_wall_tot, time_wall_tot/3600.
  write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
  &  time_cpu_tot/time_wall_tot*100.
endif
start_time_cpu = finish_time_cpu
start_time_wall = finish_time_wall
!
!
! check presence of optical components from input
if (input%xs%dfoffdiag) then
   offdiag = .true.
else
   offdiag = .false.
endif
nlf = .true.
if (input%xs%gqmax .gt. 0.d0) nlf = .false.
!
!  === do loops over all available phonon modes and displacements
!
eq_done = .false.
!
! step range off equilibrium geometry
istep_lo = -int(dble(input%properties%raman%nstep-1) / 2.d0)
istep_hi = nint(dble(input%properties%raman%nstep-1) / 2.d0)
i_shift = -istep_lo + 1
!if (input%properties%raman%nstep .le. 2) i_shift = i_shift + 1
!
! loop over optical modes
!
do imode = 1, nmode
   if ((input%properties%raman%mode .ne. 0) .and. (imode .ne. input%properties%raman%mode)) cycle
   if (rank .eq. 0) then
    write(99, '(/,"  +--------------+",/,   &
               &  "  |  Mode ",i3," :  |",/,   &
               &  "  +--------------+",/)') imode
   endif
! check for acoustic modes
   acoustic(imode) = .false.
   active(imode) = .true.
   if (input%properties%raman%usesym) then
      call check_acoustic(ev(:, imode), acoustic(imode))
      if (acoustic(imode)) then
         if (rank .eq. 0) write(99, '(/,"Info(Raman): This mode seems to be acoustic",//)')
         cycle
      endif
! check mode if Raman active by symmetry
      if (rank .eq. 0) call check_raman (imode, dble(ev(:, imode)), irep(imode), active(imode))
#ifdef MPI
      call MPI_Bcast(active, nmode, MPI_Logical, 0, MPI_Comm_World, ierr)
#endif
      if (.not. active(imode)) then
         if (rank .eq. 0) write(99, '(/,"Info(Raman): This mode is not Raman active",//)')
         cycle
      endif
#ifdef MPI
!     call MPI_Bcast(sym_rt2, 81, MPI_Double_Complex, 0, MPI_Comm_World, ierr)
#endif
   endif
!
! displace atoms along eigenvector
   i = 1
   do istep = istep_lo, istep_hi
      if (rank .eq. 0) write (99, '(/," *** Working on step ",i2," ***",/)') istep + i_shift
!
! do equilibrium geometry only once
      if (istep .eq. 0 .and. input%properties%raman%doequilibrium) then
         if (rank .eq. 0) write(99, '("Info(Raman): This is the equilibrium geometry!")')
! create unambiguous file extension and directory
         Write (raman_filext, '("_MOD", I3.3, "_DISP", I2.2, ".OUT")') 0, 0
         Write (raman_stepdir, '("./MOD", I3.3, "_DISP", I2.2, "/")') 0, 0
      else
         Write (raman_filext, '("_MOD", I3.3, "_DISP", I2.2, ".OUT")') imode, istep + i_shift
         Write (raman_stepdir, '("./MOD", I3.3, "_DISP", I2.2, "/")') imode, istep + i_shift
      endif
! check if computation was already done
      Inquire (file='RAMAN_POT'//trim(raman_filext), Exist=existent)
      if (existent) then
        if (rank .eq. 0) then
         write(99, '("Info(Raman): Groundstate calculation for mode ",i3," and step ",i2, &
           &         " seems to be already done.")') imode, istep + i_shift
         write(99, '("             Reading from file ",a)') 'RAMAN_POT'//trim(filext)
         call flushifc(99)
        endif
        call raman_readpot(read_i, 'RAMAN_POT'//trim(raman_filext), read_dph, read_engy, &
               & force_sum)
        if (read_i .ne. istep + i_shift) then
           write(*,'("Error(Raman): Step number in file ",a," not consistent!")') 'RAMAN_POT'//trim(raman_filext)
           write(*,'("Read     :  ",i2)') read_i
           write(*,'("Expected :  ",i2)') istep + i_shift
           stop
        endif
        t1 = read_dph - input%properties%raman%displ*dble(istep)
        if (istep .eq. 0 .and. .not.input%properties%raman%doequilibrium) &
    &         t1 = read_dph - 100.d0*dble(natmtot)*input%structure%epslat
        if (abs(t1) .gt. eps ) then
           write(*,'("Error(Raman): Step length in file ",a," not consistent!")') 'RAMAN_POT'//trim(raman_filext)
           write(*,'("Read     :  ",f14.6)') read_dph
           write(*,'("Expected :  ",f14.6)') input%properties%raman%displ*dble(istep)
           stop
        endif
        potinx(imode, istep + i_shift) = read_dph
        if (input%properties%raman%useforces) then
          potiny(imode, istep + i_shift) = force_sum
        else
          potiny(imode, istep + i_shift) = read_engy
        endif
      else
        if (rank .eq. 0) then
         write(99, '("Info(Raman): Performing groundstate calculation for mode ",i3," and step ",i2)') imode, i
!#ifdef IFORT
         inquire(file=trim(adjustl(raman_stepdir))//'.test', exist=existent)
         if (.not. existent) then
            j = system('mkdir '//trim(adjustl(raman_stepdir)))
            if (j .ne. 0) &
          &    write(*, '("Warning(Raman): When executing mkdir ",a,", the error ",i4," was returned. ")') &
          &         trim(adjustl(raman_stepdir)), j
            open(unit=13, file=trim(adjustl(raman_stepdir))//'.test', status='unknown')
            close(13)
         endif
!#else
!        call execute_command_line('mkdir '//trim(adjustl(raman_stepdir)))
!#endif
        endif
#ifdef MPI
        call MPI_Barrier(MPI_Comm_World, ierr)
#endif
!#ifdef IFORT
        j = chdir('./'//trim(adjustl(raman_stepdir)))
!#else
!       call execute_command_line('cd ./'//trim(adjustl(raman_stepdir)))
!#endif
! start anew from input geometry, construct distorted geometry and run through loop
        call raman_restore_struct
        call init0
        dph = input%properties%raman%displ*dble(istep)
        if (istep .eq. 0 .and. .not.input%properties%raman%doequilibrium) &
    &         dph = 100.d0*dble(natmtot)*input%structure%epslat
        call dcell(vgamc(:), ev(:, imode), dph)
! set appropriate parameters, do groundstate calculation and compute forces
        input%structure%primcell = .False.
        input%structure%autormt = .False.
        input%structure%tshift = .false.
        ! if requested in input.xml compute forces
        input%groundstate%tforce = input%properties%raman%useforces
        task = 0
        notelns = 5
        notes(1) = '                                     '
        notes(2) = '   +--------------------------------+'
        notes(3) = '   | GNDSTATE was called from RAMAN |'
        notes(4) = '   +--------------------------------+'
        notes(5) = '                                     '
        call gndstate
        if (rank .eq. 0) then
         write(99, '("             Groundstate calculation done.")')
         call flushifc(99)
        endif
        force_sum = 0.d0
        if (input%properties%raman%useforces) then
!          ! store forces times eigenvector
           do is = 1, nspecies
              do ia = 1, natoms(is)
                 ias = idxas (ia, is)
                 if (all(abs(dble(ev(:, imode))) .lt. eps)) then
                 ! account for vectors rotated to imaginary axis
                 force_sum = force_sum + (forcetot(1, ias)*aimag(ev(3*ias-2, imode)) + &
                    &                     forcetot(2, ias)*aimag(ev(3*ias-1, imode)) + &
                    &                     forcetot(3, ias)*aimag(ev(3*ias, imode)) ) / &
                    &                     sqrt(spmass(is))
                else
                 ! the usual case
                 force_sum = force_sum + (forcetot(1, ias)*dble(ev(3*ias-2, imode)) + &
                    &                     forcetot(2, ias)*dble(ev(3*ias-1, imode)) + &
                    &                     forcetot(3, ias)*dble(ev(3*ias, imode)) ) / &
                    &                     sqrt(spmass(is))
                endif
             enddo
           EndDo
           force_sum = force_sum*sqrt(fgew)
        endif
! store potential
        if (rank .eq. 0) then
         call raman_writepot(istep+i_shift, '../RAMAN_POT'//trim(raman_filext), dph, engytot, force_sum)
        endif
        potinx(imode, istep + i_shift) = dph
        if (input%properties%raman%useforces) then
          potiny(imode, istep + i_shift) = force_sum
        else
          potiny(imode, istep + i_shift) = engytot
        endif
      ! take time
        call cpu_time(t_cpu_proc)
        call timesec(finish_time_wall)
#ifdef MPI
        call MPI_Allreduce(t_cpu_proc, finish_time_cpu,  1, MPI_Real8, MPI_Sum, MPI_Comm_World, ierr)
#else
        finish_time_cpu = t_cpu_proc
#endif
        time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
        time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
        if (rank .eq. 0) then
         write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') &
           &  time_cpu_tot, time_cpu_tot/3600.
         write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') &
           &  time_wall_tot, time_wall_tot/3600.
         write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
           &  time_cpu_tot/time_wall_tot*100.
        endif
        start_time_cpu = finish_time_cpu
        start_time_wall = finish_time_wall
! clean-up
!       if (.not. associated(input%gw)) call raman_delgndst
!#ifdef IFORT
        j = chdir('..')
!#else
!       call execute_command_line('cd ..')
!#endif
      endif
!
! --------- do XS calculation ---------------
!
! check if computation was already done
      existent = .true.
      do oct1 = 1, 3
         do oct2 = 1, 3
            if (oct1 .ne. oct2 .and. .not. offdiag) cycle
            if (input%xs%xstype .eq. 'TDDFT') then
               Call genfilname (basename='EPSILON', asc=.False., &
                & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
                & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=nlf, &
                & fxctypestr=input%xs%tddft%fxctype, tq0=.true., &
                & oc1=oct1, oc2=oct2, iqmt=1, filnam=fneps)
            else
              fcoup = input%xs%bse%coupling
              epsilondir='EPSILON'
              if(fcoup) then
                tdastring=''
              else
                if(input%xs%bse%chibarq) then 
                  tdastring="-TDA-BAR"
                else
                  tdastring="-TDA"
                end if
              end if
              if(input%xs%bse%bsetype == "IP") then
                tdastring=''
              end if
              bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
              scrtypestring = '-'//trim(input%xs%screening%screentype)

               Call genfilname (dirname=trim(epsilondir),basename='EPSILON', tq0=.True., oc1=oct1, &
                & oc2=oct2, bsetype=trim(bsetypestring), &
                & scrtype=trim(scrtypestring), nar= .Not. &
                & input%xs%bse%aresbse, filnam=fneps)
            endif
            fileeps(oct1, oct2) = fneps
            inquire (file=trim(raman_stepdir)//fileeps(oct1, oct2), Exist=existent1)
            existent = existent .and. existent1
         enddo
      enddo
      if (existent) then
      ! read previously saved data
         if (rank .eq. 0) then
          write (99, '("Info(Raman): XS calculation for mode ",i3," and step ",i2,&
            &          " seems to be already done.")') imode, istep + i_shift
         endif
         do oct1 = 1, 3
            do oct2 = 1, 3
               if (oct1 .ne. oct2 .and. .not. offdiag) cycle
               if (rank .eq. 0) then
                write (99, '("             Reading dielectric function from file ",a)') &
                 &          trim(adjustl(raman_stepdir))//trim(adjustl(fileeps(oct1, oct2)))
                call flushifc(99)
               endif
               call raman_readeps (imode, istep+i_shift, oct1, oct2, trim(raman_stepdir)//fileeps(oct1, oct2))     
               if (input%properties%raman%molecule) then
               ! use molecular polarizability instead
               ! df contains (electronic) polarizability of the molecule in Bohr^3 molecule^-1
                  if (oct1 .eq. oct2) then
                     df(imode, istep+i_shift, oct1, oct2, :) = omega*(df(imode, istep+i_shift, oct1, oct2, :)-zone)/4.d0/pi
                  else
                     df(imode, istep+i_shift, oct1, oct2, :) = omega*df(imode, istep+i_shift, oct1, oct2, :)/4.d0/pi
                  endif
               endif
            enddo
         enddo
      else
         if (rank .eq. 0) &
         & write (99, '("Info(Raman): Performing XS calculation for mode ",i3," and step ",i2)') imode, istep + i_shift
#ifdef MPI
         call MPI_Barrier(MPI_Comm_World, ierr)
#endif
!#ifdef IFORT
         j = chdir('./'//trim(raman_stepdir))
!#else
!        call execute_command_line('cd ./'//trim(raman_stepdir))
!#endif
!
! launch xs part according to element xs in input.xml:
!
! start anew from input geometry, construct distorted geometry and run xs
         call raman_restore_struct
         call init0
         dph = input%properties%raman%displ*dble(istep)
         if (istep .eq. 0 .and. .not.input%properties%raman%doequilibrium) &
    &         dph = 100.d0*dble(natmtot)*input%structure%epslat
         call dcell(vgamc(:), ev(:, imode), dph)
         input%structure%primcell = .False.
         input%structure%autormt = .False.
         input%structure%tshift = .false.
         notelns = 5
         notes(1) = '                                     '
         notes(2) = '   +--------------------------------+'
         notes(3) = '   |    XS was called from RAMAN    |'
         notes(4) = '   +--------------------------------+'
         notes(5) = '                                     '
         If (associated(input%gw)) then
            if (rank .eq. 0) then
             write (99, '("Info(Raman): Performing GW calculation...")')
             call flushifc(99)
            endif
            Call gwtasklauncher
            if (rank .eq. 0) then
             write (99, '("             GW done")')
            endif
         ! take time
            call cpu_time(t_cpu_proc)
#ifdef MPI
            call MPI_Allreduce(t_cpu_proc, finish_time_cpu, 1, MPI_Real8, MPI_Sum, MPI_Comm_World, ierr)
#else
            finish_time_cpu = t_cpu_proc
#endif
            call timesec(finish_time_wall)
            time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
            time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
            if (rank .eq. 0) then
             write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') &
               &   time_cpu_tot, time_cpu_tot/3600.
             write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') &
               &   time_wall_tot, time_wall_tot/3600.
             write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
               &   time_cpu_tot/time_wall_tot*100.
             call flushifc(99)
            endif
            start_time_cpu = finish_time_cpu
            start_time_wall = finish_time_wall
         endif
!
!  compute optical spectra according to input
         call xstasklauncher
!
!  if the user did not specify a plan element, deallocate the plan constructed in xstasklauncher       
         if (.not. dplan) deallocate(input%xs%plan)
!  read from files
         do oct1 = 1, 3
            do oct2 = 1, 3
               if (oct1 .ne. oct2 .and. .not. offdiag) cycle
               call raman_readeps (imode, istep+i_shift, oct1, oct2, trim(fileeps(oct1, oct2)))
               if (input%properties%raman%molecule) then
               ! use molecular polarizability instead
               ! df contains (electronic) polarizability of the molecule in Bohr^3 molecule^-1
                  if (oct1 .eq. oct2) then
                     df(imode, istep+i_shift, oct1, oct2, :) = &
                      & omega*(df(imode, istep+i_shift, oct1, oct2, :)-zone)/4.d0/pi
                  else
                     df(imode, istep+i_shift, oct1, oct2, :) = &
                      & omega*df(imode, istep+i_shift, oct1, oct2, :)/4.d0/pi
                  endif
               endif
            enddo
         enddo
!#ifdef IFORT
         j = chdir('..')
!#else
!         call execute_command_line('cd ..')
!#endif
      endif
      if (rank .eq. 0) write(99, '("             XS calculation done.")')
      ! take time
      call cpu_time(t_cpu_proc)
#ifdef MPI
      call MPI_Allreduce(t_cpu_proc, finish_time_cpu, 1, MPI_Real8, MPI_Sum, MPI_Comm_World, ierr)
#else
      finish_time_cpu = t_cpu_proc
#endif
      call timesec(finish_time_wall)
      time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
      time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
      if (rank .eq. 0) then
       write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') &
           &  time_cpu_tot, time_cpu_tot/3600.
       write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') &
           &  time_wall_tot, time_wall_tot/3600.
       write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
           &  time_cpu_tot/time_wall_tot*100.
       call flushifc(99)
      endif
      start_time_cpu = finish_time_cpu
      start_time_wall = finish_time_wall
! turn off computation of equilibrium geometry after first mode
      if (istep .eq. 0) eq_done = .true.
      i = i + 1
   enddo
   if (rank .eq. 0) then
    write(99, '("Info(Raman): Done with mode ",i3,/)') imode
    call flushifc(99)
   endif
enddo
!
call raman_restore_struct
!
if (rank .eq. 0) write(99, '(/,"Info(Raman): *** Done with all modes! ***")')
!
! === end loops over all phonon modes and displacements
!
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ START RAMAN +++ START RAMAN +++ START RAMAN +++
! +++++++++++++++++++++++++++++++++++++++++++++++++++
!
! ok, all necessary data are hopefully collected, let's do the Raman calculation now...
!
!
if (rank .eq. 0) then
 write(99, '(/,"Info(Raman): Starting calculation of Raman scattering intensities...",/,&
          &    "             See file RAMAN.OUT for results.",/)')
 call flushifc(99)
else
 ! actual Raman calculation is not particularly heavy and hence skipped by child procs
 goto 477
endif
! number of knots for computations of potential, one interval contains 4 knots
maxp = 4*input%properties%raman%ninter + 1
!
! open main output file
open(unit=66,file='RAMAN.OUT',status='unknown',form='formatted')
!
! === start loop over all phonon modes
!
do imode = 1, nmode
   if (acoustic(imode)) cycle
   if ((input%properties%raman%mode .ne. 0) .and. (imode .ne. input%properties%raman%mode)) cycle
   if (.not. active(imode)) cycle
!
   write(66,'(116("*"),/,40("*"),"   RAMAN INTENSITIES FOR MODE ",i3,"   ",40("*"),/,116("*"),//)') imode
!
!  change file extension
   Write (filext, '("_MOD", I3.3, ".OUT")') imode
!
   open(unit=77,file='RAMAN_POTENTIAL'//trim(filext),status='unknown',form='formatted')
!
!  first fit desired functions to given data points
!  for the potential
   call polyfit(imode)
!  ...and the dielectric function
   call polyfit_diel (imode, rlas)
!
!  write PARAMETERS to OUTPUT file
   write(66,'(//,116("*"),/46("*"),"   START CALCULATION   ",47("*"))')
   write(66,'(/" Potential  coefficients:    ",7f11.5)') a0,a1,a2,a3,a4,a5,a6
   write(66,'(/" x between ",f7.4," and ",f7.4,"  Bohr")') input%properties%raman%xmin, &
             &                                             input%properties%raman%xmax
   write(66,'(/," Number of unit cells: ",2x,i8,5x," Volume of unit cell [ Bohr^3 ]: ",f16.8)')  ncell, omega
   write(66,'(  "                                                             [ cm^3 ]: ",g16.8)') omega*fau3cm3
   write(66,'(" Laser energy [ cm-1 ]: ",10x,f12.2)') fhawn*rlas
   write(66,'(/," Derivatives of the dielectric function deps/du: ",/, &
    &           " Re                                         Im")')
   write(66,'(/,"  ( ",3f12.3," )       ( ",3f12.3," ) ",/, &
    &           "  ( ",3f12.3," )       ( ",3f12.3," ) ",/, &
    &           "  ( ",3f12.3," )       ( ",3f12.3," ) ")') &
    &   (dble(deq(1, oct2)),oct2=1,3), (aimag(deq(1, oct2)),oct2=1,3), &
    &   (dble(deq(2, oct2)),oct2=1,3), (aimag(deq(2, oct2)),oct2=1,3), &
    &   (dble(deq(3, oct2)),oct2=1,3), (aimag(deq(3, oct2)),oct2=1,3)
!  endif
!
   write(66, '(/," Include local field effects for dielectric function : ",l1)') .not.nlf
   write(66, '(/," Broadening [ cm-1 ] : ",f7.2)') input%properties%raman%broad
!
   call getfgew ( ev(:, imode) )
   sfact = 0.5d0 / fgew
!  write(66, '(/," The effective mass of this mode is [ amu ]            : ",f12.3)') fgew/famuau
!
   sn = dsqrt(dble(ncell))
!
!
!    find MINIMUM of potential, version 3 using lapack
   zmin = 0.d0
   call findmin3(zmin)
!
!    SHIFT everything into MINIMUM
   zs = zmin
   call shift(zs)
   input%properties%raman%xmin = input%properties%raman%xmin - zs
   input%properties%raman%xmax = input%properties%raman%xmax - zs
   call epsshift(zs)
   write(66,'(/,37("*"),"  Phonon potential shifted into minimum  ",38("*"),/)')
   write(66,'(/" Potential  coefficients:    ",7f11.5)') a0,a1,a2,a3,a4,a5,a6
   write(66,'(/" x between ",f7.4," and ",f7.4,"  Bohr")') input%properties%raman%xmin, &
    &                                                      input%properties%raman%xmax
!
   write(66,'(/," Derivatives of the dielectric function deps/du: ",/, &
    &           " Re                                         Im")')
   write(66,'(/,"  ( ",3f12.3," )       ( ",3f12.3," ) ",/, &
    &           "  ( ",3f12.3," )       ( ",3f12.3," ) ",/, &
    &           "  ( ",3f12.3," )       ( ",3f12.3," ) ")') &
    &   (dble(deq(1, oct2)),oct2=1,3), (aimag(deq(1, oct2)),oct2=1,3), &
    &   (dble(deq(2, oct2)),oct2=1,3), (aimag(deq(2, oct2)),oct2=1,3), &
    &   (dble(deq(3, oct2)),oct2=1,3), (aimag(deq(3, oct2)),oct2=1,3)
!
!
!    determine coefficients for N cells
   if( .not. input%properties%raman%molecule) call ncells(sn,ncell)
!
!    calculate SPECTRA
!
   write(77,*) '# effective potential'
   call potential(maxp)
   call eigenen
   call transmat(input%properties%raman%ninter,h) 
   if (input%properties%raman%molecule) then
      call spectrum_m(rlas, filext)
      write(66,'(/,"Info(Raman): Spectrum computed.",/)')
   else
      do oct1 = 1, 3
         do oct2 = oct1, 3
            if (oct1 .ne. oct2 .and. .not. offdiag) cycle
            call spectrum(rlas, oct1, oct2, comp(oct1,oct2), filext)
            write(66,'(/,"Info(Raman): Spectrum for optical component ", &
  &                   a2," computed.",/)') comp(oct1,oct2)
         enddo
      enddo
   endif
   close (77)
! resonance behavior of the fundamental transition
   ! conversion |deps/du|^2 -> Raman susceptibility |dchi'/dQ|^2
   if (input%properties%raman%molecule) then
      t1 = 1.d0/omega/fgew
   else
      t1 = omega/fgew/(4.d0*pi)**2
   endif
   do oct1 = 1, 3
      do oct2 = oct1, 3
         if (oct1 .ne. oct2 .and. .not. offdiag) cycle
         open(unit=77,file='RAMAN_RESONANCE_OC'//comp(oct1,oct2)//trim(filext),&
        &             status='unknown',form='formatted')
         dwlas = (input%xs%energywindow%intv(2) - &
   &              input%xs%energywindow%intv(1)) / &
   &              dble(input%xs%energywindow%points)
         if (input%properties%raman%molecule) then
            write(77,'("# Resonance behavior ",/, &
       &          "# Elas [ Ha ]     Elas [ eV ]     Elas [ nm ]     Esc [ Ha ]     ", &
       &          "   |dchi/dQ|^2 [ au ]    ", &
       &          "(d sigma)/(d Omega) [ 10^-36 m^2 sr^-1 ]      ", &
       &          "(d sigma)/(d Omega)/Esc^4 [ 10^-60 m^6 sr^-1 ]",/)')
         else
            write(77,'("# Resonance behavior ",/, &
       &          "# Elas [ Ha ]     Elas [ eV ]     Elas [ nm ]     Esc [ Ha ]     ", &
       &          "   |dchi/dQ|^2 [ au ]    ", &
       &          "S_tot [ 10^-5 sr^-1 m^-1 ]      S_tot/Esc^4 [ 10^-29 m^3 sr^-1 ]",/)')
         endif
         do iw = 1, input%xs%energywindow%points
            wlas = input%xs%energywindow%intv(1) + dble(iw - 1)*dwlas
            if (wlas .lt. eps) cycle
            call polyfit_diel_res(imode, iw)
            call epsshift(zs)
            call spectrum_res(wlas, oct1, oct2, Sab, ws)
            write(77,'(f12.5,3x,f12.3,3x,f12.2,3x,f12.5,5x,3g20.8)') &
       &    wlas, wlas*fhaev, 1.d0/(wlas*fharnm), ws, t1*abs(deq(oct1, oct2))**2, &
       &       Sab, Sab/(ws*fhawn)**3/(wlas*fhawn)*1.d16
         enddo
         close(77)
      enddo
   enddo
!
! end loop over modes
enddo
!
!
! add CPU time used for execution
call cpu_time(finish_time_cpu)
call timesec(finish_time_wall)
write(66,'(//,116("*"),//," CPU time used in RAMAN: ",1x,f12.2," seconds",/, &
              &           " wall time used in RAMAN: ",f12.2," seconds",/)') &
              &       finish_time_cpu - start_time_cpu, finish_time_wall - start_time_wall
write(66,'(/,53("*"),"   END   ",54("*"),/116("*"))')
!
write(99,'("Info(Raman): Raman calculation done",/)')
write(99,'(" CPU time used in RAMAN: ",f12.2," seconds (",f6.2," hrs)",/)')&
&  finish_time_cpu - start_time_cpu, (finish_time_cpu - start_time_cpu)/3600.
! write end
time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') time_cpu_tot, time_cpu_tot/3600.
write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') time_wall_tot, time_wall_tot/3600.
write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
&  time_cpu_tot/time_wall_tot*100.
!
! shut down Raman specific files and arrays
close(66); close(80); close(99)
!
! child procs jump here
477 continue
!
! turn off phonons after Raman
if (associated(input%phonons)) nullify (input%phonons)
!
! turn off GW calculations after Raman
if (associated(input%gw)) nullify (input%gw)
!
! turn off XS calculations after Raman
nullify (input%xs)
!
deallocate( df )
deallocate( eigen,T,R,z1,z2 )
deallocate( transme1,transme2,transme3,transme4,transme5,transme6 )
deallocate( e1,e2,e3,de,indexi )
deallocate( xa,xpot,pot,b0,b1,b2,b3,b4 )
deallocate( potinx,potiny )
deallocate( active, acoustic, irep)
if (allocated(dynq)) deallocate(dynq)
if (allocated(dynp)) deallocate(dynp)
if (allocated(dynr)) deallocate(dynr)
!
!
!
!
end subroutine raman
!
