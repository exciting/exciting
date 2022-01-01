!
!
!
! Copyright (C) 2021 J. Cimurs, A. Gulans
!
!
Subroutine calcbfb (ikp, vnlvv)
      Use modmain
      Use modinput
      Use modgw, only : kqset,Gkqset, kset, nomax, numin, ikvbm, ikcbm, ikvcm, Gset
      Use modmpi, only : rank
      Use modfvsystem
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (In) :: vnlvv (nstsv, nstsv)

! local variables
      Integer :: ngknr, ik, jk, ist1, ist2, ist3, is4
      Integer :: is, ia, ias, ic, m1, m2, lmax, lm, ilm, irc ! ilm and lm potentially redundant (one of them)
      Integer :: nrc, iq, ig, iv (3), igq0, igk
      Integer :: ilo, loindex
      Integer :: info
      Integer :: nmatp

      Real (8) :: v (3), cfq, ta,tb
      Complex (8) zrho01, zrho02, zt1, zt2
      Integer :: nr, l, m, io1, lm2, ir, if3
      type(evsystem) :: system
! automatic arrays
      Real (8) :: zn (nspecies)
      Complex (8) sfacgq0 (natmtot)
! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: zm (:, :)
      Complex (8), Allocatable :: matrixl(:, :)
      Complex (8), Allocatable :: matrixm(:, :)
      Complex (8), Allocatable :: matrixm1(:, :)
      Complex (8), Allocatable :: matrixm2(:, :)
      Complex (8), Allocatable :: matrixPC(:, :)
      type (WFType) :: wf1,wf2,prod,pot
! external functions
      Complex (8) zfinp, zfmtinp
      External zfinp, zfmtinp

! -- Adaptively Compressed Exchange Operator starts --
      Allocate (matrixl(nstsv,nstsv))

      if (.not.(allocated(pace))) then
        Allocate (pace(nstsv, nmatmax, nkpt))
        pace=zzero
      endif
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))

      allocate(evecfv(nmat(1,ikp),nstfv))
      allocate(zm(nstfv,nmat(1,ikp)))

      nmatp=nmat(1,ikp)

      matrixl = -vnlvv ! create copy

! Remove upper triangular part
      Do ist1 = 2, nstsv
        matrixl(ist1-1,ist1:)=0
      End Do

      Call zpotrf('L',nstsv,matrixl,nstsv,info) ! Computes the Cholesky factorization
      If (info.ne.0) Then
        Write (*, *) 'vnlvv is not negative definite! Info=',info
        stop
      End If

      Call match (ngk(1, ikp), gkc(:, 1, ikp), tpgkc(:, :, 1, ikp), sfacgk(:, :, 1, ikp), apwalm)

      call newsystem(system,input%groundstate%solver%packedmatrixstorage,nmatp)
      call MTRedirect(mt_hscf%main,mt_hscf%spinless)
      call overlapsetup(system, ngk(1, ikp), apwalm, igkig(:, 1, ikp), vgkc(:,:,1,ikp))

      ! c
      call getevecfv(vkl(:,ikp), vgkl(:,:,:,ikp), evecfv)

      call zgemm('c', 'n', nstfv, nmatp, nmatp, &
      &          zone, evecfv, nmatp, &
      &          system%overlap%za, nmatp, &
      &          zzero, zm, nstfv)

      ! Vnl*conjg(c)*S
      call zgemm('n', 'n', nstfv, nmatp, nstfv, &
      &          zone, matrixl, nstfv, &
      &          zm, nstfv, zzero, &
      &          pace(1,1,ikp), nstfv)

      call deletesystem(system)
      deallocate(apwalm,zm,evecfv)

!-- Adaptively Compressed Exchange Operator ends --

      Return
End Subroutine
!EOC
