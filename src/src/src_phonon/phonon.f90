
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine phonon
      Use modmain
      use modmpi
      Use modinput
#ifdef IFORT
      use ifport
#endif
      Implicit None
! local variables
      Integer :: is, js, ia, ja, ka, jas, kas
      Integer :: iq, ip, jp, nph, iph, i
      integer :: task0
      Real (8) :: dph, a, b, t1
      Real (8) :: ftp (3, maxatoms, maxspecies)
      real(8) :: ftp0(3, maxatoms, maxspecies)
      Complex (8) zt1, zt2
      Complex (8) dyn (3, maxatoms, maxspecies)
      character(256) :: status
      logical :: finished
      logical :: gammap
      Real (8) :: genrad (3)
! allocatable arrays
      Real (8), Allocatable :: veffmtp (:, :, :)
      Real (8), Allocatable :: veffirp (:)
      Complex (8), Allocatable :: dveffmt (:, :, :)
      Complex (8), Allocatable :: dveffir (:)
! subdirectory handling
#ifndef IFORT
      integer system, chdir
#endif
      integer :: j
      logical :: existent
!------------------------!
!     initialisation     !
!------------------------!
! require forces
      input%groundstate%tforce = .True.
! no primitive cell determination
      input%structure%primcell = .False.
! use autokpt=.true. for supercells, if not requested in input file
! compute radkpt from input ngridk in this case
      if (.not. input%groundstate%autokpt) then
         genrad(:) = dble(input%groundstate%ngridk(:) - 1) * &
                     Sqrt(input%structure%crystal%basevect(1, :)**2 + &
                          input%structure%crystal%basevect(2, :)**2 + &
                          input%structure%crystal%basevect(3, :)**2)
         input%groundstate%radkpt = Max(genrad(1), genrad(2), genrad(3))
         input%groundstate%autokpt = .True.
      endif
! initialise universal variables
      Call init0
! initialise q-point dependent variables
      Call init2
! write q-points to file
      Open (50, File='QPOINTS_PHONON.OUT', Action='WRITE', Form='FORMATTED')
      Write (50, '(I6, " : nqpt; q-point, vql, wqpt below")') &
     & nqpt
      Do iq = 1, nqpt
         Write (50, '(I6, 4G18.10, 2I8)') iq, vql (:, iq), wqpt (iq)
      End Do
      Close (50)
! read original effective potential from file and store in global arrays
      Call readstate
      If (allocated(veffmt0)) deallocate (veffmt0)
      Allocate (veffmt0(lmmaxvr, nrmtmax, natmtot))
      If (allocated(veffir0)) deallocate (veffir0)
      Allocate (veffir0(ngrtot))
      veffmt0 (:, :, :) = veffmt (:, :, :)
      veffir0 (:) = veffir (:)
! allocate local arrays
      Allocate (dveffmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (dveffir(ngrtot))
! switch off automatic determination of muffin-tin radii
      input%structure%autormt = .False.
! no shifting of atomic basis allowed
      input%structure%tshift = .False.
! store original parameters
      natoms0 (1:nspecies) = natoms (1:nspecies)
      natmtot0 = natmtot
      avec0 (:, :) = input%structure%crystal%basevect(:, :)
      ainv0 (:, :) = ainv (:, :)
      atposc0 (:, :, :) = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            atposc0 (:, ia, is) = atposc (:, ia, is)
         End Do
      End Do
      ngrid0 (:) = ngrid (:)
      ngrtot0 = ngrtot
!---------------------------------------!
!     compute dynamical matrix rows     !
!---------------------------------------!
      if (input%phonons%gamma .eq. 'onestep') then
        ! re-run the ground-state calculation for forces at equlibirum geo
        task0 = task
        task = 1
        Call gndstate
        ! store the total force for the equlibrium geometry
        Do js = 1, nspecies
           Do ja = 1, natoms (js)
              jas = idxas (ja, js)
              ftp0 (:, ja, js) = forcetot (:, jas)
           End Do
        End Do
        task = task0
      endif
!
10    Continue
      natoms (1:nspecies) = natoms0 (1:nspecies)
! find a dynamical matrix to calculate
      if (rank.eq.0) then
        Call dyntask (iq, is, ia, ip, status)
        finished=.false.
        if (status .eq. "finished") finished=.true.
      end if
#ifdef MPI
      Call MPI_bcast (iq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Call MPI_bcast (ia, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Call MPI_bcast (is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Call MPI_bcast (ip, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Call MPI_bcast (finished, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
      if (finished) goto 20
! phonon dry run
      If (task .Eq. 201) Go To 10
! check to see if mass is considered infinite
      If (spmass(is) .Le. 0.d0) Then
         if (rank .eq. 0) then
            Call phfext (iq, is, ia, ip, 0, 1, filext, filextdyn, phdirname)
            Open (80, File='DYN'//trim(filextdyn), Action='WRITE', Form='FORMATTED')
            Do ip = 1, 3
               Do js = 1, nspecies
                  Do ja = 1, natoms0 (js)
                     Do jp = 1, 3
                        Write (80, '(2G18.10, " : is = ", I4, ", ia = ", I4, ", ip = ", I4)') 0.d0, 0.d0, js, ja, jp
                     End Do
                  End Do
               End Do
            End Do
            Close (80)
         endif
         Go To 10
      End If
      task = 200
      nph = 1
      gammap = .false.
      If (all(ivq(:, iq) .Eq. 0)) then
         nph = 0
         gammap = .true.
      endif
      dyn (:, :, :) = 0.d0
      dveffmt (:, :, :) = 0.d0
      dveffir (:) = 0.d0
! loop over phases (cos and sin displacements)
      Do iph = 0, nph
! restore input values
         natoms (1:nspecies) = natoms0 (1:nspecies)
         input%structure%crystal%basevect(:, :) = avec0 (:, :)
         atposc (:, :, :) = atposc0 (:, :, :)
! generate the supercell
         Call phcell (iph, input%phonons%deltaph, iq, is, ia, ip)
! generate file names
         Call phfext (iq, is, ia, ip, iph, 1, filext, filextdyn, phdirname)
! create and enter subdirectory
         if (rank .eq. 0) then
            inquire(file=trim(adjustl(phdirname))//'.test', exist=existent)
            if (.not. existent) then
               j = system('mkdir '//trim(adjustl(phdirname)))
               if (j .ne. 0) &
             &    write(*, '("Warning(Phonon): When executing mkdir ",a,", the error ",i4," was returned. ")') &
             &         trim(adjustl(phdirname)), j
               open(unit=13, file=trim(adjustl(phdirname))//'.test', status='unknown')
               close(13)
            endif
         endif
#ifdef MPI
        call MPI_Barrier(MPI_Comm_World, ierr)
#endif
         j = chdir('./'//trim(adjustl(phdirname)))
! run the ground-state calculation
         if ( (rank .eq. 0) .and. (task .eq. 1) ) call writestate
         Call gndstate
! read STATE.OUT file with current extension
         call readstate
         if (rank.eq.0)  Call phdelete
! store the total force for the first displacement
         Do js = 1, nspecies
            Do ja = 1, natoms (js)
               jas = idxas (ja, js)
               ftp (:, ja, js) = forcetot (:, jas)
            End Do
         End Do
! store the effective potential for the first displacement
         Allocate (veffmtp(lmmaxvr, nrmtmax, natmtot))
         Allocate (veffirp(ngrtot))
         veffmtp (:, :, :) = veffmt (:, :, :)
         veffirp (:) = veffir (:)
! restore input values
         natoms (1:nspecies) = natoms0 (1:nspecies)
         input%structure%crystal%basevect(:, :) = avec0 (:, :)
         atposc (:, :, :) = atposc0 (:, :, :)
         if (gammap) then
            select case (input%phonons%gamma)
            case ('onestep')
              ! nothing more to do
              goto 67
            case ('twostep')
              ! displace in negative direction
              dph = -input%phonons%deltaph
            case default
              dph = input%phonons%deltaph + input%phonons%deltaph
            end select
         else
            ! generate the supercell again with twice the displacement
            dph = input%phonons%deltaph + input%phonons%deltaph
         endif
         Call phcell (iph, dph, iq, is, ia, ip)
! generate new file names
         Call phfext (iq, is, ia, ip, iph, 2, filext, filextdyn, phdirname)
! write STATE.OUT file with updated extension for use in gndstate
         if ( rank .eq. 0) call writestate
#ifdef MPI
         call MPI_Barrier(MPI_Comm_World, ierr)
#endif
! run the ground-state calculation again starting from the previous density
         task = 1
         Call gndstate
! read STATE.OUT file with current extension
         call readstate
         if (rank.eq.0)  Call phdelete
 67      continue
! compute the complex perturbing effective potential with implicit q-phase
         Call phdveff (iph, iq, veffmtp, veffirp, dveffmt, dveffir)
         if (input%phonons%gamma .eq. 'twostep' .and. gammap) then
            dveffmt(:, :, :) = 0.5d0 * dveffmt(:, :, :)
            dveffir(:) = 0.5d0 * dveffir(:)
         endif
         Deallocate (veffmtp, veffirp)
! Fourier transform the force differences to obtain the dynamical matrix
         zt1 = 1.d0 / (dble(nphcell)*input%phonons%deltaph)
! multiply by i for sin-like displacement
         If (iph .Eq. 1) zt1 = zt1 * zi
         kas = 0
         Do js = 1, nspecies
            ka = 0
            Do ja = 1, natoms0 (js)
               Do i = 1, nphcell
                  ka = ka + 1
                  kas = kas + 1
                  t1 = - dot_product (vqc(:, iq), vphcell(:, i))
                  zt2 = zt1 * cmplx (Cos(t1), Sin(t1), 8)
                  Do jp = 1, 3
                    if (gammap) then
                      select case (input%phonons%gamma)
                      case ('onestep')
                        t1 = ftp0(jp, ka, js) - ftp(jp, ka, js)
                      case ('twostep')
                        t1 = 0.5d0 * (forcetot(jp, kas)-ftp(jp, ka, js))
                      case default
                        t1 = - (forcetot(jp, kas)-ftp(jp, ka, js))
                      end select
                    else
                      t1 = - (forcetot(jp, kas)-ftp(jp, ka, js))
                    endif
                    dyn (jp, ja, js) = dyn (jp, ja, js) + zt2 * t1
                  End Do
               End Do
            End Do
         End Do
! return to base directory
         j = chdir('..')
      End Do
! write dynamical matrix row to file
      Open (80, File='DYN'//trim(filextdyn), Action='WRITE', Form='FORMATTED')
      Do js = 1, nspecies
         Do ja = 1, natoms0 (js)
            Do jp = 1, 3
               a = dble (dyn(jp, ja, js))
               b = aimag (dyn(jp, ja, js))
               If (Abs(a) .Lt. 1.d-12) a = 0.d0
               If (Abs(b) .Lt. 1.d-12) b = 0.d0
               if (rank.eq.0) Write (80, '(2G18.10, " : is = ", I4, ", ia = ", I4, ", &
              &ip = ", I4)') a, b, js, ja, jp
            End Do
         End Do
      End Do
      Close (80)
! write the complex perturbing effective potential to file
      if (rank.eq.0)  Call writedveff (iq, is, ia, ip, dveffmt, dveffir)
      Go To 10
20 continue
End Subroutine
