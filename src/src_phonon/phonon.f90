
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine phonon
      Use modmain
      use modmpi
      Use modinput
      Implicit None
! local variables
      Integer :: is, js, ia, ja, ka, jas, kas
      Integer :: iq, ip, jp, nph, iph, i
      Real (8) :: dph, a, b, t1
      Real (8) :: ftp (3, maxatoms, maxspecies)
      Complex (8) zt1, zt2
      Complex (8) dyn (3, maxatoms, maxspecies)
      character(256) :: status
      logical :: finished
! allocatable arrays
      Real (8), Allocatable :: veffmtp (:, :, :)
      Real (8), Allocatable :: veffirp (:)
      Complex (8), Allocatable :: dveffmt (:, :, :)
      Complex (8), Allocatable :: dveffir (:)
!------------------------!
!     initialisation     !
!------------------------!
! require forces
      input%groundstate%tforce = .True.
! no primitive cell determination
      input%structure%primcell = .False.
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
! determine k-point grid size from radkpt
      ! let the user decide whether to use an automatic k-point grid or not
      !input%groundstate%autokpt = .True.
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
10    Continue
      natoms (1:nspecies) = natoms0 (1:nspecies)
! find a dynamical matrix to calculate
      if (rank.eq.0) then
        Call dyntask (80, iq, is, ia, ip, status)
        finished=.false.
        if (status .eq. "finished") finished=.true.
      end if
#ifdef MPI
      Call MPI_bcast (iq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Call MPI_bcast (ia, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Call MPI_bcast (is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Call MPI_bcast (ip, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Call MPI_bcast (finished, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      Call phfext (iq, is, ia, ip, filext)
#endif
      if (finished) goto 20
! phonon dry run
      If (task .Eq. 201) Go To 10
! check to see if mass is considered infinite
      If (spmass(is) .Le. 0.d0) Then
         Do ip = 1, 3
            Do js = 1, nspecies
               Do ja = 1, natoms0 (js)
                  Do jp = 1, 3
                     if (rank.eq.0) Write (80, '(2G18.10, " : is = ", I4, ", ia = ", I&
                    &4, ", ip = ", I4)') 0.d0, 0.d0, js, ja, jp
                  End Do
               End Do
            End Do
         End Do
         if (rank.eq.0) Close (80)
         Go To 10
      End If
      task = 200
      nph = 1
      If ((ivq(1, iq) .Eq. 0) .And. (ivq(2, iq) .Eq. 0) .And. (ivq(3, &
     & iq) .Eq. 0)) nph = 0
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
! run the ground-state calculation
         Call gndstate
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
! generate the supercell again with twice the displacement
         dph = input%phonons%deltaph + input%phonons%deltaph
         Call phcell (iph, dph, iq, is, ia, ip)
! run the ground-state calculation again starting from the previous density
         task = 1
         Call gndstate
! compute the complex perturbing effective potential with implicit q-phase
         Call phdveff (iph, iq, veffmtp, veffirp, dveffmt, dveffir)
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
                     t1 = - (forcetot(jp, kas)-ftp(jp, ka, js))
                     dyn (jp, ja, js) = dyn (jp, ja, js) + zt2 * t1
                  End Do
               End Do
            End Do
         End Do
      End Do
! write dynamical matrix row to file
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
! delete the non-essential files
      if (rank.eq.0)  Call phdelete
      Go To 10
20 continue
End Subroutine
