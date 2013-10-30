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
use mod_atoms, only: natmtot, nspecies, natoms, spmass
use mod_energy, only: engytot
use modinput
use modxs
use raman_input
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
!
implicit none
integer :: maxp,it,ntp
integer :: i, j, ia, is, iat, ic, imode, istep, nmode
integer :: istep_lo, istep_hi, i_shift
integer :: oct, oct1, oct2
integer :: read_i
real(8) :: rlas,sn,temp,zmin,zs, norm
real(8) :: dph, vgamc(3)
real(8) :: start_time_cpu, finish_time_cpu, time_cpu_tot
real(8) :: start_time_wall, finish_time_wall, time_wall_tot
real(8) :: read_dph, read_engy
Real(8), Allocatable :: w(:)
Complex(8), Allocatable :: ev(:, :)
Complex(8), Allocatable :: dynq(:, :, :)
Complex(8), Allocatable :: dynp(:, :)
Complex(8), Allocatable :: dynr(:, :, :)
Logical, Allocatable :: active(:), acoustic(:)
integer, allocatable :: irep(:)
Logical :: existent, existent1, eq_done, nlf, lt2(3, 3)
Character(256) :: raman_filext, raman_stepdir
character(80) :: ext
!
! take time
time_cpu_tot = 0.d0; time_wall_tot = 0.d0
call cpu_time(start_time_cpu)
call timesec(start_time_wall)
!
call init0
call init2
!
! default number of modes
nmode = 3*natmtot
!
!   unit no.,'./filename','old/unknown','FORMATTED'
!   5 		input file
!   66 		main output file
!   72 		spectrum x (opened in subroutine spectrum)
!   73 		spectrum y "
!   74 		spectrum z "
!   77 		potential for plotting
!   80 		dielectric function for plotting
!   250ff       if requested in input, several files containing eigenvalues and functions for plotting
!               (opened in subroutine eigenen)
!
!
! open INFO file for all pre-Raman computations triggered in this subroutine
open(unit=99,file='INFO_RAMAN.OUT',status='unknown',form='formatted')
!
!     boundaries
xmin_r = input%properties%raman%xmin
xmax_r = input%properties%raman%xmax
!
numpot = input%properties%raman%nstep
if (input%properties%raman%nstep .le. 2) then
! minimalistic determination of coefficients, equilibirum + 1 distortion gives quadratic potential and linear epsilon dep.
   input%properties%raman%nstep = 2
   numpot = 3
endif
allocate( potinx(3*natmtot, numpot) )
allocate( potiny(3*natmtot, numpot) )
!
! scattering volume for molecules or solids
if (input%properties%raman%molecule) then
   ncell = 1
else
   ncell = 100000
endif
! cell volume in cm3 (for later use with cm-1)
!     vol_cm3 = omega*1.d-24*faua**3
!
!     laser energy   
select case(trim(input%properties%raman%elaserunit))
  case ('eV')
     rlas = input%properties%raman%elaser * fevha
  case ('nm')
     rlas = 1.d0/input%properties%raman%elaser * frnmha
  case ('cm-1')
     rlas = input%properties%raman%elaser * fwnha
  case default
end select
!     rlas is now in Ha
!
!
!     broadening of fundamental and overtones
gamma1 = input%properties%raman%broad
gamma2 = input%properties%raman%broad
gamma3 = input%properties%raman%broad
gamma4 = input%properties%raman%broad
!
!
!     temperature range, lattice expansion
!     steps of tempi from tempa to tempe
tempa = input%properties%raman%temp
tempe = input%properties%raman%temp
tempi = 0.d0
!
!     read(5,*) intphonon                              ! in a solid, take other N-1 phonons into account (1) or not (0)
!     intphonon = 0
!     if (input%properties%raman%intphonon) intphonon = 1
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
! eigen values and vectors, dynamical matrices
Allocate (w(3*natmtot))
Allocate (ev(3*natmtot, 3*natmtot))
Allocate (active(3*natmtot))
Allocate (acoustic(3*natmtot))
Allocate (irep(3*natmtot))
Allocate (dynq(3*natmtot, 3*natmtot, nqpt))
Allocate (dynp(3*natmtot, 3*natmtot))
Allocate (dynr(3*natmtot, 3*natmtot, ngridq(1)*ngridq(2)*ngridq(3)))
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
!
! +++ SYMMETRY ANALYSIS OF THE CRYSTAL +++
! +++ CONSTRUCTION OF CHARACTER TABLE  +++
! 
input%structure%tshift = .true.
!call init0
call construct_chartabl
write(99,'("Info(Raman): Character table constructed.",/)')
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
      write(99,'("Info(Raman): Performing a supercell phonon calculation for the Gamma point.",/)')
      call flushifc(99)
      ! create dummy phonon element if not present in the input file
      If ( .Not. (associated(input%phonons))) Then
         input%phonons => getstructphonons (emptynode)
      End If
      If (associated(input%phonons%parts)) nullify (input%phonons%parts)
      ! set parameters for Gamma phonons
      input%phonons%do = 'fromscratch'
      input%phonons%ngridq = (/ 1, 1, 1/)
      call phonon
      ! restore original parameters
      call rereadinput
      ! continue with reading from files...
      Call readdyn (.true.,dynq)
      ! apply the acoustic sum rule
      Call sumrule (dynq)
      ! obtain eigenvectors at Gamma
      Call dynqtor (dynq, dynr)
      Call dynrtoq (vgamc(:), dynr, dynp)
      Call dyndiag (dynp, w, ev)
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
   Case ('fromfile')
   ! read in the dynamical matrix from files
      write(99,'("Info(Raman): Reading dynamical matrix from files.",/)')
      ! create dummy phonon element if not present in the input file
      If ( .Not. (associated(input%phonons))) Then
         input%phonons => getstructphonons (emptynode)
      End If
      ! initialize general and q-dependent variables
      call init2
      ! read dynamical matrix
      Call readdyn (.true.,dynq)
      ! apply the acoustic sum rule
      Call sumrule (dynq)
      ! obtain eigenvectors at Gamma
      Call dynqtor (dynq, dynr)
      Call dynrtoq (vgamc(:), dynr, dynp)
      Call dyndiag (dynp, w, ev)
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
   Case ('symvec')
   ! construct symmetry vectors and use them as normal coordinates
      call construct_symvec(ev)
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
   Case ('symveccheck')
   ! construct symmetry vectors and return
      call construct_symvec(ev)
      write(99, '(/," Symmetry vectors ",/, &
           &        " ---------------- ",/)')
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
      if (associated(input%phonons)) nullify (input%phonons)
      if (associated(input%gw)) nullify (input%gw)
      nullify (input%xs)
      return
   Case ('readinput')
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
         write(99, '("Info(Raman): Norm of eigenvector is ",f14.9)') sqrt(norm)
         write(99, '("             Eigenvector will be renormalized")')
         ev = ev / sqrt(norm)
      endif
      call flushifc(99)
      ! use input variable mode for direct loop over modes
      nmode = 1
End Select
!
! take time
call cpu_time(finish_time_cpu)
call timesec(finish_time_wall)
time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') time_cpu_tot, time_cpu_tot/3600.
write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') time_wall_tot, time_wall_tot/3600.
write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
&  time_cpu_tot/time_wall_tot*100.
start_time_cpu = finish_time_cpu
start_time_wall = finish_time_wall

!
! check presence of optical components from input
if (input%xs%xstype .eq. 'TDDFT' .and. input%xs%dfoffdiag) then
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
if (input%properties%raman%nstep .le. 2) i_shift = i_shift + 1
!
! loop over optical modes
!do imode = 4,3*natmtot
!
do imode = 1, nmode
   if ((input%properties%raman%mode .ne. 0) .and. (imode .ne. input%properties%raman%mode)) cycle
   write(99, '(/,"  +--------------+",/,   &
              &  "  |  Mode ",i3," :  |",/,   &
              &  "  +--------------+",/)') imode
   write(*,*) 'mode ',ev(:, imode)
! check for acoustic modes
   call check_acoustic(dble(ev(:, imode)), acoustic(imode))
   if (acoustic(imode)) cycle
! check mode if Raman active by symmetry
   call check_raman (dble(ev(:, imode)), irep(imode), active(imode))
   if (.not. active(imode)) then
      write(99, '(/,"Info(Raman): This mode is not Raman active",//)')
      cycle
   endif
!
!  analyze symmetry of Raman tensor
   Write (ext, '("_MOD", I3.3, ".OUT")') imode
   call construct_t2 (irep(imode), ext, lt2)
   write(*, *) lt2(1, :)
   write(*, *) lt2(2, :)
   write(*, *) lt2(3, :)
!
! displace atoms along eigenvector
   i = 1
   do istep = istep_lo, istep_hi
      write (*,  '(/," *** Working on step ",i2," ***",/)') istep + i_shift
      write (99, '(/," *** Working on step ",i2," ***",/)') istep + i_shift
!
! do equilibrium geometry only once
      if (istep .eq. 0) then
         write(99, '("Info(Raman): This is the equlibirum geometry!")')
! create unambiguous file extension and directory
         Write (raman_filext, '("_MOD", I3.3, "_DISP", I2.2, ".OUT")') 0, 0
         Write (raman_stepdir, '("./MOD", I3.3, "_DISP", I2.2, "/")') 0, 0
!        Write (scrpath, '("./MOD", I3.3, "_DISP", I2.2, "/")') 0, 0
      else
         Write (raman_filext, '("_MOD", I3.3, "_DISP", I2.2, ".OUT")') imode, istep + i_shift
         Write (raman_stepdir, '("./MOD", I3.3, "_DISP", I2.2, "/")') imode, istep + i_shift
!        Write (scrpath, '("./MOD", I3.3, "_DISP", I2.2, "/")') imode, istep + i_shift
      endif
! check if computation was already done
      Inquire (file='RAMAN_POT'//trim(raman_filext), Exist=existent)
      if (existent) then
         write(99, '("Info(Raman): Groundstate calculation for mode ",i3," and step ",i2, &
           &         " seems to be already done.")') imode, istep + i_shift
         write(99, '("             Reading from file ",a)') 'RAMAN_POT'//trim(filext)
         call raman_readpot(read_i, 'RAMAN_POT'//trim(raman_filext), read_dph, read_engy)
         if (read_i .ne. istep + i_shift) then
            write(*,'("Error(Raman): Step number in file ",a," not consistent!")') 'RAMAN_POT'//trim(raman_filext)
            write(*,'("Read     :  ",i2)') read_i
            write(*,'("Expected :  ",i2)') istep + i_shift
            stop
         endif
         potinx(imode, istep + i_shift) = read_dph
         potiny(imode, istep + i_shift) = read_engy
         call flushifc(99)
      else
         write(99, '("Info(Raman): Performing groundstate calculation for mode ",i3," and step ",i2)') imode, i
! intermediate solution using system and chdir functions
!#ifdef IFORT
         j = system('mkdir '//trim(raman_stepdir))
         j = chdir('./'//trim(raman_stepdir))
         j = system('cp ../input.xml .')
!#else
!         call execute_command_line('mkdir '//trim(raman_stepdir))
!         call execute_command_line('cd ./'//trim(raman_stepdir))
!         call execute_command_line('cp ../input.xml .')
!#endif
! start anew from input geometry, construct distorted geometry and run through loop
         call rereadinput
         call init0
         write(*, '("Equilibrium geometry coordinates")')
         do is = 1, nspecies
            Do ia = 1, natoms (is)
               write(*,*) 'is, ia, coord(:) ', is, ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            enddo
         enddo
         dph = input%properties%raman%displ*dble(istep)
         write(*,'("Phonon eigenvector")')
         iat = 1
         do is = 1, nspecies
            Do ia = 1, natoms (is)
               write(*,*) 'is, ia, eigvec(:) ', is, ia, dble(ev((3*(iat-1)+1):(3*iat), imode))
               iat = iat + 1
            enddo
         enddo
         call dcell(vgamc(:), ev(:, imode), dph)
         write(*, '("             Supercell constructed:")' )
         do is = 1, nspecies
            Do ia = 1, natoms (is)
               write(*,*) 'is, ia, coord(:) ', is, ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            enddo
         enddo
!        call writegeometryraman
! set appropriate parameters, do groundstate calculation and compute forces
         input%structure%primcell = .False.
         input%structure%autormt = .False.
         input%structure%tshift = .false.
         input%groundstate%tforce = .true.
         task = 0
         notelns = 5
         notes(1) = '                                     '
         notes(2) = '   +--------------------------------+'
         notes(3) = '   | GNDSTATE was called from RAMAN |'
         notes(4) = '   +--------------------------------+'
         notes(5) = '                                     '
         write(*,*) 'call to groundstate'
         call flushifc(99)
         call gndstate
         write(*,*) 'groundstate computed'
         write(99, '("             Groundstate calculation done.")')
         call flushifc(99)
! store forces
!        Do jas = 1, natmtot
!           force_sum = force_sum + forcetot (:, jas)
!        End Do
! store potential
         write(*,*) dph, engytot
         call raman_writepot(istep+i_shift, '../RAMAN_POT'//trim(raman_filext), dph, engytot)
         potinx(imode, istep + i_shift) = dph
         potiny(imode, istep + i_shift) = engytot
      ! take time
         call cpu_time(finish_time_cpu)
         call timesec(finish_time_wall)
         time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
         time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
         write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') &
           &  time_cpu_tot, time_cpu_tot/3600.
         write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') &
           &  time_wall_tot, time_wall_tot/3600.
         write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
           &  time_cpu_tot/time_wall_tot*100.
         start_time_cpu = finish_time_cpu
         start_time_wall = finish_time_wall
! clean-up
!        if (.not. associated(input%gw)) call raman_delgndst
!#ifdef IFORT
         j = chdir('..')
!#else
!         call execute_command_line('cd ..')
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
            Inquire (file='EPSILON_OC'//comp(oct1,oct2)//trim(raman_filext), Exist=existent1)
            existent = existent .and. existent1
         enddo
      enddo
      if (existent) then
      ! read previously saved data
         write (99, '("Info(Raman): XS calculation for mode ",i3," and step ",i2,&
           &          " seems to be already done.")') imode, istep + i_shift
         do oct1 = 1, 3
            do oct2 = 1, 3
               if (oct1 .ne. oct2 .and. .not. offdiag) cycle
               write (99, '("             Reading dielectric function from file ",a)')  &
                &         'EPSILON_OC'//comp(oct1,oct2)//trim(filext)
               call raman_readeps (imode, istep+i_shift, oct1, oct2, 'EPSILON_OC'//comp(oct1,oct2)//trim(raman_filext))
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
         call flushifc(99)
      else
         write (99, '("Info(Raman): Performing XS calculation for mode ",i3," and step ",i2)') imode, istep + i_shift
!#ifdef IFORT
         j = chdir('./'//trim(raman_stepdir))
         j = system('cp ../input.xml .')
!#else
!         call execute_command_line('cd ./'//trim(raman_stepdir))
!         call execute_command_line('cp ../input.xml .')
!#endif
!
! launch xs part according to element xs in input.xml:
!
! start anew from input geometry, construct distorted geometry and run xs
         call rereadinput
         call init0
         dph = input%properties%raman%displ*dble(istep)
         call dcell(vgamc(:), ev(:, imode), dph)
         write(*, '("             Supercell constructed:")' )
         do is = 1, nspecies
            Do ia = 1, natoms (is)
               write(*,*) 'is, ia, coord(:) ', is, ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            enddo
         enddo
!        call writegeometryraman
! set appropriate parameters, do groundstate calculation and compute forces
         input%structure%primcell = .False.
         input%structure%autormt = .False.
         input%structure%tshift = .false.
         notelns = 5
         notes(1) = '                                     '
         notes(2) = '   +--------------------------------+'
         notes(3) = '   |    XS was called from RAMAN    |'
         notes(4) = '   +--------------------------------+'
         notes(5) = '                                     '
         call flushifc(99)
         If (associated(input%gw)) then
            write (99, '("Info(Raman): Performing GW calculation...")')
            call flushifc(99)
            Call gwtasklauncher
            write (99, '("             GW done")')
         ! take time
            call cpu_time(finish_time_cpu)
            call timesec(finish_time_wall)
            time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
            time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
            write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') &
              &   time_cpu_tot, time_cpu_tot/3600.
            write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') &
              &   time_wall_tot, time_wall_tot/3600.
            write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
              &   time_cpu_tot/time_wall_tot*100.
            start_time_cpu = finish_time_cpu
            start_time_wall = finish_time_wall
            call flushifc(99)
         endif
         write(*,*) ' call to XS'
         call xstasklauncher
         write(*,*) 'XS computed'
         if (input%xs%xstype .eq. 'TDDFT') then
            do oct1 = 1, 3
               Do oct2 = 1, 3
                  if (oct1 .ne. oct2 .and. .not. offdiag) cycle
                  Call genfilname (basename='EPSILON', asc=.False., &
                 & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
                 & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=nlf, &
                 & fxctypestr=input%xs%tddft%fxctype, tq0=.true., &
                 & oc1=oct1, oc2=oct2, iqmt=1, filnam=fneps)
                  call raman_readeps (imode, istep + i_shift, oct1, oct2, fneps)
                  call raman_writeeps (imode, istep + i_shift, oct1, oct2, comp(oct1,oct2), raman_filext)
                  if (input%properties%raman%molecule) then
                  ! use molecular polarizability instead
                  ! df contains (electronic) polarizability of the molecule in Bohr^3 molecule^-1
                     if (oct1 .eq. oct2) then
                        df(imode, istep+i_shift, oct1, oct2, :) = omega*(df(imode, istep+i_shift, oct1, oct2, :)-zone)/4.d0/pi
                     else
                        df(imode, istep+i_shift, oct1, oct2, :) = omega*df(imode, istep+i_shift, oct1, oct2, :)/4.d0/pi
                     endif
                  endif
               Enddo
            Enddo
         else ! it is BSE
            Do oct = 1, 3
               Call genfilname (basename='EPSILON', tq0=.True., oc1=oct, &
              & oc2=oct, bsetype=input%xs%bse%bsetype, &
              & scrtype=input%xs%screening%screentype, nar= .Not. &
              & input%xs%bse%aresbse, filnam=fneps)
               call raman_readeps (imode, istep + i_shift, oct, oct, fneps)
               call raman_writeeps (imode, istep + i_shift, oct, oct, comp(oct,oct), raman_filext)
               if (input%properties%raman%molecule) then
               ! use molecular polarizability instead
               ! df contains (electronic) polarizability of the molecule in Bohr^3 molecule^-1
                  df(imode, istep+i_shift, oct, oct, :) = omega*(df(imode, istep+i_shift, oct, oct, :)-zone)/4.d0/pi
               endif
            Enddo
         endif
! copy relevant files
!        call copyxsfiles
! clean-up
!        call raman_delxs
!#ifdef IFORT
         j = chdir('..')
!#else
!         call execute_command_line('cd ..')
!#endif
      endif
      write(99, '("             XS calculation done.")')
      ! take time
      call cpu_time(finish_time_cpu)
      call timesec(finish_time_wall)
      time_cpu_tot = time_cpu_tot + finish_time_cpu - start_time_cpu
      time_wall_tot = time_wall_tot + finish_time_wall - start_time_wall
      write(99,'(/," Total CPU time used: ",1x,f12.2," seconds (",f9.2," hrs)")') &
          &  time_cpu_tot, time_cpu_tot/3600.
      write(99,'(" Total wall time used: ",f12.2," seconds (",f9.2," hrs)")') &
          &  time_wall_tot, time_wall_tot/3600.
      write(99,'(" Av. CPU utilization: ",f12.2," percent ",/)')&
          &  time_cpu_tot/time_wall_tot*100.
      start_time_cpu = finish_time_cpu
      start_time_wall = finish_time_wall
      call flushifc(99)
! turn off computation of equilibrium geometry after first mode
      if (istep .eq. 0) eq_done = .true.
      i = i + 1
   enddo
! mirror data for minimalistic nstep=2 case
   if (input%properties%raman%nstep .le. 2) then
      potinx(imode, 1) = potinx(imode, 3)
      potiny(imode, 1) = potiny(imode, 3)
      df(imode, 1, :, :, :) = ztwo*df(imode, 2, :, :, :) - df(imode, 3, :, :, :)
   endif
   write(99, '("Info(Raman): Done with mode ",i3,/)') imode
   call flushifc(99)
enddo
write(99, '(/,"Info(Raman): *** Done with all modes! ***")')
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
write(99, '(/,"Info(Raman): Starting calculation of Raman scattering intensities...",/,&
         &    "             See file RAMAN.OUT for results.",/)')
call flushifc(99)
 
maxp = 4*input%properties%raman%ninter + 1           ! number of knots for computations of potential, one interval contains 4 knots
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
      call polyfit (imode)
      write(*, '("Info(Raman): Potential fitted")')
!  ...and the dielectric function
      call polyfit_diel (imode, rlas)
      write(*, '("Info(Raman): Dielectric function fitted")')
!
!  write PARAMETERS to OUTPUT file
      write(66,'(//,116("*"),/46("*"),"   START CALCULATION   ",47("*"))')
      write(66,'(/" Potential  coefficients:    ",7f11.5)') a0,a1,a2,a3,a4,a5,a6
      write(66,'(/" x between ",f7.4," and ",f7.4,"  Bohr")') xmin_r,xmax_r
      write(66,'(/," Number of unit cells: ",2x,i8,5x," Volume of unit cell [ Bohr^3 ]: ",f16.8)')  ncell, omega
      write(66,'(  "                                                             [ cm^3 ]: ",g16.8)') omega*fau3cm3
      write(66,'(" Laser energy [ cm-1 ]: ",10x,f12.2)') fhawn*rlas
      do oct1 = 1, 3
         Do oct2 = 1, 3
            if (oct1 .ne. oct2 .and. .not. offdiag) cycle
            write(66,'(/," Derivatives of the dielectric function for optical component ",a2,/,&
     &                   "Der",11x,"Re",15x,"Im",/,6(i3,f17.3,f17.3,/))') &
     &                   comp(oct1, oct2),1, deq(oct1, oct2),&
     &                                    2,d2eq(oct1, oct2),&
     &                                    3,d3eq(oct1, oct2),&
     &                                    4,d4eq(oct1, oct2),&
     &                                    5,d5eq(oct1, oct2),&
                                          6,d6eq(oct1, oct2)
         enddo
      enddo
      write(66, '(/," Include local field effects for dielectric function : ",l1)') .not.nlf
      write(66, '(/," Broadening [ cm-1 ] : ",4f7.2)') gamma1,gamma2,gamma3,gamma4
!
      call getfgew ( ev(:, imode) )
      sfact = 0.5d0 / fgew
!
      if (abs(tempe - tempa) .lt. 1.d-5) then
         ntp = 1
      else
         if (abs(tempi) .lt. 1.d-5) then
            ntp = 1
         else
            ntp = int( (tempe - tempa)/tempi ) + 1                ! number of temp steps
         endif
      endif
!     if ((intphonon .eq. 1) .and. (ntp .gt. 1)) then
!        write(*,*) 'Warning! Phonon interaction for temp loops currently not available'
!        write(66,*) ' Phonon interaction for temperature loop switched off.'
!        intphonon = 0 
!     endif
      sn = dsqrt(dble(ncell))
!     num1 = nnumber
!     if (num2 .eq. 0) num2 = 1
!
!
!
!  interaction with other phonons, if requested
!     if (intphonon .eq. 1) then
!        call potential(maxp)
!        call eigenen(num1,.false.)
!        write(*,*) ' preliminary eigen solver finished'
!
!      determine COEFFICIENTS for TEMPERATURE temp
!        call temppot(temp)
!        call potential(maxp)
!        write(*,*) ' temp averaged multi-phonon potential computed'
!     endif
!
!    find MINIMUM of potential, version 3 using lapack
      zmin = 0.d0
      call findmin3(zmin)
!
!    SHIFT everything into MINIMUM
      zs = zmin
      call shift(zs)
      xmin_r = xmin_r - zs
      xmax_r = xmax_r - zs
      call epsshift(zs)
      write(*,*) ' shift to potential minimum done'
      write(66,'(/,37("*"),"  Phonon potential shifted into minimum  ",38("*"),/)')
      write(66,'(/" Potential  coefficients:    ",7f11.5)') a0,a1,a2,a3,a4,a5,a6
      write(66,'(/" x between ",f7.4," and ",f7.4,"  Bohr")') xmin_r,xmax_r
      do oct1 = 1, 3
         Do oct2 = 1, 3
            if (oct1 .ne. oct2 .and. .not. offdiag) cycle
            write(66,'(/," Derivatives of the dielectric function for optical component ",a2,/,&
     &                   "Der",11x,"Re",15x,"Im",/,6(i3,f17.3,f17.3,/))') &
     &                   comp(oct1, oct2),1, deq(oct1, oct2),&
     &                                    2,d2eq(oct1, oct2),&
     &                                    3,d3eq(oct1, oct2),&
     &                                    4,d4eq(oct1, oct2),&
     &                                    5,d5eq(oct1, oct2),&
                                          6,d6eq(oct1, oct2)
         enddo
      enddo

!
!    determine coefficients for N cells
      if( .not. input%properties%raman%molecule) call ncells(sn,ncell)
!
!    calculate SPECTRA
!
      write(77,*) '# effective potential'
      call potential(maxp)
      call eigenen
      write(*,*) ' main eigen solver finished'
      call transmat(input%properties%raman%ninter,h) 
      write(*,*) ' vib matrix elements computed'
!    TEMPERATURE loop
      do it = 1, ntp
         temp = tempa + (it - 1)*tempi
         if (input%properties%raman%molecule) then
            call spectrum_m(temp, rlas, filext)
            write(66,'(/,"Info(Raman): Spectrum computed.",/)')
         else
            do oct1 = 1, 3
               do oct2 = oct1, 3
                  if (oct1 .ne. oct2 .and. .not. offdiag) cycle
                  call spectrum(temp, rlas, oct1, oct2, comp(oct1,oct2), filext)
                  write(66,'(/,"Info(Raman): Spectrum for optical component ", &
     &                      a2," computed.",/)') comp(oct1,oct2)
               enddo
            enddo
         endif
      enddo
      close (77)
!
! end loop over modes
enddo
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
!
!
!
!
end subroutine raman
!
