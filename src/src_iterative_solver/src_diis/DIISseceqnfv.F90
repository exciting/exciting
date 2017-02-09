
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine DIISseceqnfv (ik, ispn, apwalm, vgpc, evalfv, evecfv)
!
  !USES:
      Use modmain, Only: nstfv, vkl, ngk, igkig, nmat, vgkl, timemat, &
     & npmat, apwordmax, lmmaxapw, natmtot, nkpt, nmatmax, nspnfv, &
     & timefv, ngkmax, zzero, zone
      Use sclcontroll
      Use diisinterfaces
      Use modfvsystem
  ! !INPUT/OUTPUT PARAMETERS:
  !   ik     : k-point number (in,integer)
  !   ispn   : first-variational spin index (in,integer)
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   vgpc   : G+k-vectors in Cartesian coordinates
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  ! This routine will perform several Bock Davidson iterations following the sceme:
  ! 1. for each of m bands:
  !   a. calculate Residual
  !$$
  !\ket{\mathbf{R}\left(\ket{\mathbf{A}^{ap}},E^{ap}\right)}=(\mathbf{H}-E^{ap}\mathbf{S})\ket{ \mathbf{A}^{ap}}
  !$$
  !   b. calculate $\delta \mathbf{A}$
  ! 2. solve Projected system in evecsv+$\delta \mathbf{A}$ subspace
  !EOP
  !BOC
      Implicit None
  ! argumentstrialvec
      Integer, Intent (In) :: ik
      Integer, Intent (In) :: ispn
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (Inout) :: evalfv (nstfv, nspnfv)
      Complex (8), Intent (Inout) :: evecfv (nmatmax, nstfv, nspnfv)
!
  ! local variables
!
      Type (evsystem) :: system
      Logical :: packed, jacdav
      Integer :: is, ia, idiis, n, np, ievec, i, info, flag, icurrent
      Real (8) :: vl, vu, abstol
      Real (8) :: cpu0, cpu1
      Real (8) :: eps, rnorm
      Complex (8), Allocatable :: P (:, :)
      Complex (8), Allocatable :: h (:, :, :)
      Complex (8), Allocatable :: s (:, :, :)
      Complex (8), Allocatable :: r (:, :)
      Complex (8), Allocatable :: trialvecs (:, :, :)
      Complex (8), Allocatable :: eigenvector (:, :)
      Real (8), Allocatable :: eigenvalue (:, :)
!
      Real (8) :: w (nmatmax), rnorms (nstfv)
      Complex (8) :: z
      Integer :: iunconverged, evecmap (nstfv)
!
      If ((ik .Lt. 1) .Or. (ik .Gt. nkpt)) Then
         Write (*,*)
         Write (*, '("Error(seceqnfv): k-point out of range : ", I8)') &
        & ik
         Write (*,*)
         Stop
      End If
      n = nmat (ik, ispn)
      np = npmat (ik, ispn)
      Allocate (P(nmatmax, nmatmax))
      Allocate (h(nmat(ik, ispn), nstfv, maxdiisspace))
      Allocate (s(nmat(ik, ispn), nstfv, maxdiisspace))
      Allocate (r(nmat(ik, ispn), nstfv))
      Allocate (trialvecs(nmat(ik, ispn), nstfv, maxdiisspace))
      Allocate (eigenvector(nmat(ik, ispn), nstfv))
      Allocate (eigenvalue(nstfv, maxdiisspace+1))
!
  !----------------------------------------!
  !     Hamiltonian and overlap set up     !
  !----------------------------------------!
      Call cpu_time (cpu0)
      packed = .False.
      jacdav = .False.
      Call newsystem (system, packed, n)
      Call hamiltonandoverlapsetup (system, ngk(ik, ispn), apwalm, &
     & igkig(1, ik, ispn), vgpc)
!
      Call cpu_time (cpu1)
!
  !$OMP CRITICAL
      timemat = timemat + cpu1 - cpu0
  !$OMP END CRITICAL
  !update eigenvectors with iteration
      recalculate_preconditioner = .False.
      Call cpu_time (cpu0)
      If (calculate_preconditioner()) Then
         P = 0
         w = 0
         Call seceqfvprecond (n, system, P, w, evalfv(:, ispn), &
        & evecfv(:, :, ispn))
         Call writeprecond (ik, n, P, w)
      Else
     !---------------------------------!
     ! initialisation from file        !
     !---------------------------------!
         iunconverged = nstfv
         Call readprecond (ik, n, P, w)
     !    write(*,*)"readeigenvalues",w
         Call getevecfv (vkl(1, ik), vgkl(1, 1, ik, 1), evecfv)
         Call getevalfv (vkl(1, ik), evalfv)
!
         Call zlarnv (2, iseed, n*nstfv, eigenvector)
         eigenvector = cmplx (dble(eigenvector), 0.)
         Call zscal (n*nstfv, dcmplx(1e-3/n/nstfv, 0.), eigenvector, 1)
     !coppy eigenvectors to work aray eigenvector
         Do i = 1, nstfv
            Call zcopy (n, evecfv(1, i, ispn), 1, eigenvector(1, i), 1)
        !     call zaxpy(n ,zone,evecfv(1,i,ispn),1,eigenvector(1,i),1)
            eigenvalue (i, 1) = evalfv (i, ispn)
            evecmap (i) = i
         End Do
!
     !initialisation for jacobidavidson preconditioning
         If (jacdav) Call jacdavblock (n, iunconverged, system, n, &
        & eigenvector, h(:, :, idiis), s(:, :, idiis), eigenvalue(:, &
        & idiis), trialvecs(:, :, idiis), h(:, :, idiis), 0)
!
!
     !#####################
     ! start diis iteration
     !#####################
!
         Do idiis = 1, diismax
            icurrent = Mod (idiis-1, maxdiisspace) + 1
            Write (*,*) "icurrent", icurrent
            Write (*,*) "diisiter", idiis
        !----------------------------------------------------!
        ! h(:,:,diis) holds matrix with current aproximate   !
        ! vectors multiplied with hamilton                   !
        ! o: same for overlap*evecfv                         !
        !----------------------------------------------------!
!
            If (idiis .Gt. 1) Then
           !after first iteration copy refined vectors to evecfv
               Do i = 1, nstfv
                  If (evecmap(i) .Ne. 0) Call zcopy (n, eigenvector(1, &
                 & evecmap(i)), 1, evecfv(1, i, ispn), 1)
               End Do
           ! setuphs computes H*v and S*v and normalises
               Call setuphsvect (n, nstfv, system, evecfv, nmatmax, &
              & h(:, :, icurrent), s(:, :, icurrent))
           !orthogonalise all vectors
               Call orthogonalise (n, nstfv, evecfv(:, :, ispn), &
              & nmatmax, s(:, :, icurrent))
           !copy back orthogonalised vectors to work array
               Do i = 1, nstfv
                  If (evecmap(i) .Ne. 0) Call zcopy (n, evecfv(1, i, &
                 & ispn), 1, eigenvector(1, evecmap(i)), 1)
               End Do
            End If
        ! setuphs computes H*v and S*v and normalises
            Call setuphsvect (n, iunconverged, system, eigenvector, n, &
           & h(:, :, icurrent), s(:, :, icurrent))
            Call rayleighqotient (n, iunconverged, eigenvector, h(:, :, &
           & icurrent), s(:, :, icurrent), eigenvalue(:, icurrent))
            Call residualvectors (n, iunconverged, h(:, :, icurrent), &
           & s(:, :, icurrent), eigenvalue(:, icurrent), r, rnorms)
        !update eigenvalues
            Do i = 1, nstfv
               If (evecmap(i) .Ne. 0) evalfv (i, ispn) = eigenvalue &
              & (evecmap(i), icurrent)
            End Do
        !------------------------------------------------------------------!
        !check for convergence and remove converged vectors from iteration !
        !------------------------------------------------------------------!
            If (allconverged(iunconverged, rnorms) .Or. idiis .Eq. &
           & (diismax-1)) Exit
            Call remove_converged (evecmap, iunconverged, rnorms, n, r, &
           & h, s, eigenvector, eigenvalue, trialvecs)
            If (rnorms(idamax(iunconverged, rnorms, 1)) .Gt. 1e-1 .And. &
           & (idiis .Gt. 1)) Then
               recalculate_preconditioner = .True.
               Write (*,*) "recalculate preconditioner"
               Exit
           !----------------------------------------------------!
           !if all residuals are converged exit diis loop!      !
           !vectors are normalized and already copied to evecfv !
           !----------------------------------------------------!
            End If
        !-----------------------------------------------------!
        ! correction equation with spectral precond or jacdav !
        !-----------------------------------------------------!
            If ( .Not. jacdav) Then
               Call calcupdatevectors (n, iunconverged, P, w, r, &
              & eigenvalue(:, icurrent), eigenvector, trialvecs(:, :, &
              & icurrent))
               Call setuphsvect (n, iunconverged, system, trialvecs(:, &
              & :, icurrent), n, h(:, :, icurrent), s(:, :, icurrent))
               Call zcopy (n*iunconverged, trialvecs(1, 1, icurrent), &
              & 1, eigenvector, 1)
            Else
           !  call jacdavblock(n, iunconverged, system, n, &
           !  eigenvector, h(:,:,icurrent), s(:,:,icurrent), eigenvalue(:,icurrent), &
           !  trialvecs(:,:,icurrent), h(:,:,icurrent), 1)
           !  call zaxpy(n*iunconverged,zone,trialvecs(1,1,icurrent),1,eigenvector(1,1),1)
           !  call zcopy(n*iunconverged,trialvecs(1,1,icurrent),1,eigenvector(1,1),1)
            End If
        !-----------------!
        ! diis refinement !
        !-----------------!
            If (idiis .Gt. 1) Then
               Call diisupdate (idiis, icurrent, iunconverged, n, h, s, &
              & trialvecs, eigenvalue, eigenvector, info)
            End If
        !----------------!
        ! end DIIS cycle !
        !----------------!
         End Do
     !--------------------------------------!
     ! if failed recalculate preconditioner !
     !--------------------------------------!
         If (recalculate_preconditioner .Or. (idiis .Gt. diismax-1)) &
        & Then
            Call seceqfvprecond (n, system, P, w, evalfv(:, ispn), &
           & evecfv(:, :, ispn))
            Call writeprecond (ik, n, P, w)
            Write (*,*) "recalculate preconditioner"
         End If
         Call cpu_time (cpu1)
     !if(jacdav)     call jacdavblock(n, iunconverged, system, n, &
     !    eigenvector(:,idiis), h(:,:,idiis), s(:,:,idiis), eigenvalue(:,idiis), &
     !   trialvecs(:,:,idiis), h(:,:,idiis), -1)
!
      End If
!
      Call deletesystem (system)
      Deallocate (eigenvalue)
      Deallocate (eigenvector)
      Deallocate (trialvecs)
      Deallocate (r)
      Deallocate (s)
      Deallocate (h)
      Deallocate (P)
!
      timefv = timefv + cpu1 - cpu0
!
      Return
End Subroutine DIISseceqnfv
!EOC
