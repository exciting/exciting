!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: init0
! !INTERFACE:
!
!
Subroutine init0
! !USES:
      Use modinput
      Use modmain
      Use modxcifc
#ifdef XS
      Use modxs
#endif
! !DESCRIPTION:
!   Performs basic consistency checks as well as allocating and initialising
!   global variables not dependent on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, js, ia, ias
      Integer :: ist, l, m, lm, iv (3)
      Real (8) :: ts0, ts1, tv3 (3)
!
!-------------------------------!
!     zero timing variables     !
!-------------------------------!
      timeinit = 0.d0
      timemat = 0.d0
      timefv = 0.d0
      timesv = 0.d0
      timerho = 0.d0
      timepot = 0.d0
      timefor = 0.d0
      Call timesec (ts0)
!
!------------------------------------!
!     angular momentum variables     !
!------------------------------------!
      lmmaxvr = (input%groundstate%lmaxvr+1) ** 2
      lmmaxapw = (input%groundstate%lmaxapw+1) ** 2
      lmmaxmat = (input%groundstate%lmaxmat+1) ** 2
      lmmaxinr = (input%groundstate%lmaxinr+1) ** 2
      If (input%groundstate%lmaxvr .Gt. input%groundstate%lmaxapw) Then
         Write (*,*)
         Write (*, '("Error(init0): lmaxvr > lmaxapw : ", 2I8)') &
        & input%groundstate%lmaxvr, input%groundstate%lmaxapw
         Write (*,*)
         Stop
      End If
      If (input%groundstate%lmaxmat .Gt. input%groundstate%lmaxapw) &
     & Then
         Write (*,*)
         Write (*, '("Error(init0): lmaxmat > lmaxapw : ", 2I8)') &
        & input%groundstate%lmaxmat, input%groundstate%lmaxapw
         Write (*,*)
         Stop
      End If
! index to (l,m) pairs
      If (allocated(idxlm)) deallocate (idxlm)
      Allocate (idxlm(0:input%groundstate%lmaxapw,-&
     & input%groundstate%lmaxapw:input%groundstate%lmaxapw))
      lm = 0
      Do l = 0, input%groundstate%lmaxapw
         Do m = - l, l
            lm = lm + 1
            idxlm (l, m) = lm
         End Do
      End Do
! array of i**l values
      If (allocated(zil)) deallocate (zil)
      Allocate (zil(0:input%groundstate%lmaxapw))
      Do l = 0, input%groundstate%lmaxapw
         zil (l) = zi ** l
      End Do
!
!------------------------------------!
!     index to atoms and species     !
!------------------------------------!
! check if the system is an isolated molecule
      If (input%structure%molecule) Then
         input%structure%primcell = .False.
         input%structure%tshift = .False.
      End If
! find primitive cell if required
      If (input%structure%primcell) Call findprim
      natmmax = 0
      ias = 0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = ias + 1
            idxas (ia, is) = ias
         End Do
! maximum number of atoms over all species
         natmmax = Max (natmmax, natoms(is))
      End Do
! total number of atoms
      natmtot = ias
!
!------------------------!
!     spin variables     !
!------------------------!
      If (isspinspiral()) Then
         Select Case (task)
         Case (2, 3, 15, 51, 52, 53, 61, 62, 63, 120, 121)
            Write (*,*)
            Write (*, '("Error(init0): spin-spirals do not work with ta&
           &sk ", I4)') task
            Write (*,*)
            Stop
         End Select
         If (input%groundstate%xctypenumber .Lt. 0) Then
            Write (*,*)
            Write (*, '("Error(init0): spin-spirals do not work with th&
           &e OEP method")')
            Write (*,*)
            Stop
         End If
      End If
! spin-orbit coupling or fixed spin moment implies spin-polarised calculation
!
! number of spinor components and maximum allowed occupancy
      If (associated(input%groundstate%spin)) Then
         nspinor = 2
         occmax = 1.d0
      Else
         nspinor = 1
         occmax = 2.d0
      End If
! number of spin-dependent first-variational functions per state
      If (isspinspiral()) Then
         nspnfv = 2
      Else
         nspnfv = 1
      End If
! spin-polarised calculations require second-variational eigenvectors
      If (associated(input%groundstate%spin)) input%groundstate%tevecsv &
     & = .True.
! Hartree-Fock/RDMFT requires second-variational eigenvectors
      If ((task .Eq. 5) .Or. (task .Eq. 6) .Or. (task .Eq. 300)) &
     & input%groundstate%tevecsv = .True.
! get exchange-correlation functional data
      Call getxcdata ( xctype, xcdescr, xcspin, &
     & xcgrad)
      If ((associated(input%groundstate%spin)) .And. (xcspin .Eq. 0)) &
     & Then
         Write (*,*)
         Write (*, '("Error(init0): requested spin-polarised run with s&
        &pin-unpolarised")')
         Write (*, '(" exchange-correlation functional")')
         Write (*,*)
         Stop
      End If
! check for collinearity in the z-direction and set the dimension of the
! magnetisation and exchange-correlation vector fields
      If (associated(input%groundstate%spin)) Then
         ndmag = 1
         If ((Abs(input%groundstate%spin%bfieldc(1)) .Gt. &
        & input%structure%epslat) .Or. &
        & (Abs(input%groundstate%spin%bfieldc(2)) .Gt. &
        & input%structure%epslat)) ndmag = 3
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               If ((Abs(input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(1)) .Gt. input%structure%epslat) .Or. &
              & (Abs(input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(2)) .Gt. input%structure%epslat)) ndmag = 3
            End Do
         End Do
! source-free fields and spin-spirals are non-collinear in general
         If ((input%groundstate%nosource) .Or. (isspinspiral())) ndmag &
        & = 3
! spin-orbit coupling is non-collinear in general
         If (isspinorb()) ndmag = 3
      Else
         ndmag = 0
      End If
! set the non-collinear flag
      If (ndmag .Eq. 3) Then
         ncmag = .True.
      Else
         ncmag = .False.
      End If
      If ((ncmag) .And. (xcgrad .Gt. 0)) Then
         Write (*,*)
         Write (*, '("Warning(init0): GGA inconsistent with non-colline&
        &ar magnetism")')
      End If
! set fixed spin moment effective field to zero
      bfsmc (:) = 0.d0
! set muffin-tin FSM fields to zero
      bfsmcmt (:, :, :) = 0.d0
!
!-------------------------------------!
!     lattice and symmetry set up     !
!-------------------------------------!
! generate the reciprocal lattice vectors and unit cell volume
      Call reciplat
! compute the inverse of the lattice vector matrix
      Call r3minv (input%structure%crystal%basevect, ainv)
! compute the inverse of the reciprocal vector matrix
      Call r3minv (bvec, binv)
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
! map atomic lattice coordinates to [0,1) if not in molecule mode
            If ( .Not. input%structure%molecule) Call r3frac (input%structure%epslat, &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), iv)
! determine atomic Cartesian coordinates
            Call r3mv (input%structure%crystal%basevect, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), &
           & atposc(:, ia, is))
! lattice coordinates of the muffin-tin magnetic fields
            Call r3mv (ainv, input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:), bflmt(:, ia, is))
         End Do
      End Do
! lattice coordinates of the global magnetic field
      If (associated(input%groundstate%spin)) Then
         tv3 = input%groundstate%spin%bfieldc
      Else
         tv3 = 0
      End If
      Call r3mv (ainv, tv3, bfieldl)
! Cartesian coordinates of the spin-spiral vector
      If (associated(input%groundstate%spin)) Then
         tv3 = input%groundstate%spin%vqlss
      Else
         tv3 = 0
      End If
      Call r3mv (bvec, tv3, vqcss)
! find Bravais lattice symmetries
      Call findsymlat
! use only the identity if required
      If (input%groundstate%nosym) nsymlat = 1
! find the crystal symmetries and shift atomic positions if required
      Call findsymcrys
! find the site symmetries
      Call findsymsite
#ifdef XS
! determine inverse symmery elements
      Call findsymi (input%structure%epslat, maxsymcrys, nsymcrys, &
     & symlat, lsplsymc, vtlsymc, isymlat, scimap)
! generate symmetrization array for rank 2 tensors
      Call gensymt2 (maxsymcrys, nsymcrys, symlatc, lsplsymc, symt2)
! calculate advanced information on symmetry group
      Call setupsym
#endif
! automatically determine the muffin-tin radii if required
      If (input%structure%autormt) Call autoradmt
! check for overlapping muffin-tins
      Call checkmt
!
!-----------------------!
!     radial meshes     !
!-----------------------!
      nrmtmax = 1
      nrcmtmax = 1
      js = 1
      Do is = 1, nspecies
! make the muffin-tin mesh commensurate with lradstp
         nrmt (is) = nrmt (is) - Mod (nrmt(is)-1, &
        & input%groundstate%lradstep)
         nrmtmax = Max (nrmtmax, nrmt(is))
! number of coarse radial mesh points
         nrcmt (is) = (nrmt(is)-1) / input%groundstate%lradstep + 1
         nrcmtmax = Max (nrcmtmax, nrcmt(is))
! smallest muffin-tin radius
         If (rmt(is) .Lt. rmt(js)) js = is
      End Do
      If ((input%groundstate%isgkmax .Lt. 1) .Or. &
     & (input%groundstate%isgkmax .Gt. nspecies)) &
     & input%groundstate%isgkmax = js
! set up atomic and muffin-tin radial meshes
      Call genrmesh
!
!--------------------------------------!
!     charges and number of states     !
!--------------------------------------!
      chgzn = 0.d0
      chgcr = 0.d0
      chgval = 0.d0
      spnstmax = 0
      Do is = 1, nspecies
! nuclear charge
         chgzn = chgzn + spzn (is) * dble (natoms(is))
! find the maximum number of atomic states
         spnstmax = Max (spnstmax, spnst(is))
! compute the electronic charge for each species, as well as the total core and
! valence charge
         spze (is) = 0.d0
         Do ist = 1, spnst (is)
            spze (is) = spze (is) + spocc (ist, is)
            If (spcore(ist, is)) Then
               chgcr = chgcr + dble (natoms(is)) * spocc (ist, is)
            Else
               chgval = chgval + dble (natoms(is)) * spocc (ist, is)
            End If
         End Do
      End Do
! add excess charge
      chgval = chgval + input%groundstate%chgexs
! total charge
      chgtot = chgcr + chgval
      If (chgtot .Lt. 1.d-8) Then
         Write (*,*)
         Write (*, '("Error(init0): zero total charge")')
         Write (*,*)
         Stop
      End If
! effective Wigner radius
      rwigner = (3.d0/(fourpi*(chgtot/omega))) ** (1.d0/3.d0)

!-------------------------!
!     G-vector arrays     !
!-------------------------!
! find the G-vector grid sizes
      Call gridsize
! generate the G-vectors
      Call gengvec
! generate the spherical harmonics of the G-vectors
      Call genylmg
! allocate structure factor array for G-vectors
      If (allocated(sfacg)) deallocate (sfacg)
      Allocate (sfacg(ngvec, natmtot))
! generate structure factors for G-vectors
      Call gensfacgp (ngvec, vgc, ngvec, sfacg)
! generate the characteristic function
      Call gencfun
!
!-------------------------!
!     atoms and cores     !
!-------------------------!
#ifdef XS
      If (init0symonly) Go To 10
#endif
! solve the Kohn-Sham-Dirac equations for all atoms
      Call allatoms
! allocate core state eigenvalue array and set to default
      If (allocated(evalcr)) deallocate (evalcr)
      Allocate (evalcr(spnstmax, natmtot))
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist = 1, spnst (is)
               evalcr (ist, ias) = speval (ist, is)
            End Do
         End Do
      End Do
! allocate core state radial wavefunction array
      If (allocated(rwfcr)) deallocate (rwfcr)
      Allocate (rwfcr(spnrmax, 2, spnstmax, natmtot))
! allocate core state charge density array
      If (allocated(rhocr)) deallocate (rhocr)
      Allocate (rhocr(spnrmax, natmtot))
#ifdef XS
10    Continue
#endif
!
!---------------------------------------!
!     charge density and potentials     !
!---------------------------------------!
! allocate charge density arrays
      If (allocated(rhomt)) deallocate (rhomt)
      Allocate (rhomt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(rhoir)) deallocate (rhoir)
      Allocate (rhoir(ngrtot))
! allocate magnetisation arrays
      If (allocated(magmt)) deallocate (magmt)
      If (allocated(magir)) deallocate (magir)
      If (associated(input%groundstate%spin)) Then
         Allocate (magmt(lmmaxvr, nrmtmax, natmtot, ndmag))
         Allocate (magir(ngrtot, ndmag))
      End If
! Coulomb potential
      If (allocated(vclmt)) deallocate (vclmt)
      Allocate (vclmt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(vclir)) deallocate (vclir)
      Allocate (vclir(ngrtot))
! exchange-correlation potential
      If (allocated(vxcmt)) deallocate (vxcmt)
      Allocate (vxcmt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(vxcir)) deallocate (vxcir)
      Allocate (vxcir(ngrtot))
! exchange-correlation magnetic field
      If (allocated(bxcmt)) deallocate (bxcmt)
      If (allocated(bxcir)) deallocate (bxcir)
      If (associated(input%groundstate%spin)) Then
         Allocate (bxcmt(lmmaxvr, nrmtmax, natmtot, ndmag))
         Allocate (bxcir(ngrtot, ndmag))
      End If
! exchange energy density
      If (allocated(exmt)) deallocate (exmt)
      Allocate (exmt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(exir)) deallocate (exir)
      Allocate (exir(ngrtot))
! correlation energy density
      If (allocated(ecmt)) deallocate (ecmt)
      Allocate (ecmt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(ecir)) deallocate (ecir)
      Allocate (ecir(ngrtot))
! effective potential
      If (allocated(veffmt)) deallocate (veffmt)
      Allocate (veffmt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(veffir)) deallocate (veffir)
      Allocate (veffir(ngrtot))
      If (allocated(veffig)) deallocate (veffig)
      Allocate (veffig(ngvec))
! allocate muffin-tin charge and moment arrays
      If (allocated(chgmt)) deallocate (chgmt)
      Allocate (chgmt(natmtot))
      If (allocated(mommt)) deallocate (mommt)
      Allocate (mommt(3, natmtot))
!
!--------------------------------------------!
!     forces and structural optimisation     !
!--------------------------------------------!
      If (allocated(forcehf)) deallocate (forcehf)
      Allocate (forcehf(3, natmtot))
      If (allocated(forcecr)) deallocate (forcecr)
      Allocate (forcecr(3, natmtot))
      If (allocated(forceibs)) deallocate (forceibs)
      Allocate (forceibs(3, natmtot))
      If (allocated(forcetot)) deallocate (forcetot)
      Allocate (forcetot(3, natmtot))
      If (allocated(forcetp)) deallocate (forcetp)
      Allocate (forcetp(3, natmtot))
      If (allocated(tauatm)) deallocate (tauatm)
      Allocate (tauatm(natmtot))
! initialise the previous force
      forcetp (:, :) = 0.d0
! initial step sizes
      If (associated(input%structureoptimization)) Then
         tauatm (:) = input%structureoptimization%tau0atm
      Else
         tauatm (:) = 0
      End If
!
!-------------------------!
!     LDA+U variables     !
!-------------------------!
      If ((ldapu .Ne. 0) .Or. (task .Eq. 17)) Then
! LDA+U requires second-variational eigenvectors
         input%groundstate%tevecsv = .True.
! density matrices
         If (allocated(dmatlu)) deallocate (dmatlu)
         Allocate (dmatlu(lmmaxlu, lmmaxlu, nspinor, nspinor, natmtot))
! potential matrix elements
         If (allocated(vmatlu)) deallocate (vmatlu)
         Allocate (vmatlu(lmmaxlu, lmmaxlu, nspinor, nspinor, natmtot))
! zero the potential
         vmatlu (:, :, :, :, :) = 0.d0
! energy for each atom
         If (allocated(engyalu)) deallocate (engyalu)
         Allocate (engyalu(natmtot))
! interpolation constants (alpha)
         If (allocated(alphalu)) deallocate (alphalu)
         Allocate (alphalu(natmtot))
      End If
!
!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine the nuclear-nuclear energy
      Call energynn
! get smearing function data
      Call getsdata (input%groundstate%stypenumber, sdescr)
! generate the spherical harmonic transform (SHT) matrices
      Call genshtmat
!
! allocate 1D plotting arrays
      If (allocated(dvp1d)) deallocate (dvp1d)
      Allocate (dvp1d(nvp1d))
      If (allocated(vplp1d)) deallocate (vplp1d)
      Allocate (vplp1d(3, npp1d))
      If (allocated(dpp1d)) deallocate (dpp1d)
      Allocate (dpp1d(npp1d))
! zero self-consistent loop number
      iscl = 0
      tlast = .False.
!
      Call timesec (ts1)
      timeinit = timeinit + ts1 - ts0
!
      Return
End Subroutine
!EOC
