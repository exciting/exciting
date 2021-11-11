!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine oepvnlk4 (ikp, vnlcv, vnlvv)
      Use modmain
      Use modinput
      Use modgw, only : kqset,Gkqset, kset, nomax, numin, ikvbm, ikcbm, ikvcm, Gset
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (Out) :: vnlcv (ncrmax, natmtot, nstsv)
      Complex (8), Intent (Out) :: vnlvv (nstsv, nstsv)
! local variables
      Integer :: ngknr, ik, jk, ist1, ist2, ist3, is4
      Integer :: is, ia, ias, ic, m1, m2, lmax, lm
      Integer :: nrc, iq, ig, iv (3), igq0
      Real (8) :: v (3), cfq, ta,tb
      Complex (8) zrho01, zrho02, zt1, zt2
! automatic arrays
      Real (8) :: zn (nspecies)
      Complex (8) sfacgq0 (natmtot)
! allocatable arrays
      Integer, Allocatable :: igkignr (:)
      Real (8), Allocatable :: vgklnr (:, :)
      Real (8), Allocatable :: vgkcnr (:, :)
      Real (8), Allocatable :: gkcnr (:)
      Real (8), Allocatable :: tpgkcnr (:, :)
      Real (8), Allocatable :: vgqc (:, :)
      Real (8), Allocatable :: tpgqc (:, :)
      Real (8), Allocatable :: gqc (:)
      Real (8), Allocatable :: jlgqr (:, :, :)
      Real (8), Allocatable :: jlgq0r (:, :, :)
      Real (8), Allocatable :: evalsvp (:)
      Real (8), Allocatable :: evalsvnr (:)
      Real (8), Allocatable :: evalfv (:,:)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: sfacgknr (:, :)
      Complex (8), Allocatable :: ylmgq (:, :)
      Complex (8), Allocatable :: sfacgq (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :, :, :)
      Complex (8), Allocatable :: wfir1 (:, :, :)
      Complex (8), Allocatable :: wfir2 (:, :, :)
      Complex (8), Allocatable :: wfcr1 (:, :, :)
      Complex (8), Allocatable :: wfcr2 (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :, :)
      Complex (8), Allocatable :: zvclir (:, :)
      Complex (8), Allocatable :: zvcltp (:, :)
      Complex (8), Allocatable :: zfmt (:, :)

      type (WFType) :: wf1,wf2,prod,pot
! external functions
      Complex (8) zfinp, zfmtinp
      External zfinp, zfmtinp
! allocate local arrays
      if (.not.allocated(gntyyy)) call gengntyyy
      Allocate (igkignr(ngkmax))
      Allocate (vgklnr(3, ngkmax))
      Allocate (vgkcnr(3, ngkmax))
      Allocate (gkcnr(ngkmax))
      Allocate (tpgkcnr(2, ngkmax))
      Allocate (vgqc(3, ngvec))
      Allocate (tpgqc(2, ngvec))
      Allocate (gqc(ngvec))
      Allocate (jlgqr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1, ngvec, nspecies))
      Allocate (jlgq0r(0:input%groundstate%lmaxvr, nrcmtmax, nspecies))
!      Allocate (evalsvp(nstsv))
      Allocate (evalsvnr(nstsv))
      Allocate (sfacgknr(ngkmax, natmtot))
!      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
!      Allocate (evecfv(nmatmax, nstfv))
!      Allocate (evecsv(nstsv, nstsv))
      Allocate (ylmgq(lmmaxvr, ngvec))
      Allocate (sfacgq(ngvec, natmtot))
      Allocate (wfmt1(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfir1(ngrtot, nspinor, nstsv))
      Allocate (wfir2(ngrtot, nspinor, nstsv))
      Allocate (wfcr1(lmmaxvr, nrcmtmax, 2))
      Allocate (wfcr2(lmmaxvr, nrcmtmax, 2))
      Allocate (zrhomt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvclmt(lmmaxvr, nrcmtmax, natmtot, nstsv))
      Allocate (zvclir(ngrtot, nstsv))
      Allocate (zvcltp(lmmaxvr, nrcmtmax))
      Allocate (zfmt(lmmaxvr, nrcmtmax))

    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,kset%nkpt))
    evalfv(:,:) = 0.d0
    do ik = 1, nkpt
      call getevalfv(kset%vkl(:,ik), evalfv(:,ik))
    end do
    call find_vbm_cbm(1, nstfv, kset%nkpt, evalfv, efermi, nomax, numin, ikvbm, ikcbm, ikvcm)


      call WFInit(wf1)
      call WFInit(wf2)
      call WFInit(prod)
      call WFInit(pot)

!      call genWF(ikp,wf1)
!      call genWFinMT(wf1)
!      call genWFonMesh(wf1)

      allocate(pot%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
      allocate(pot%ir(ngrtot,1))


! factor for long-range term
      cfq = 0.5d0 * (omega/pi) ** 2
! set the nuclear charges to zero
      zn (:) = 0.d0
      zvclir (:, :) = 0.d0
      zvclmt (:, :, :, :) = 0.d0
      vnlcv (:, :, :) = 0.d0
      vnlvv (:, :) = 0.d0

      zrhomt (:, :, :) = 0.d0
      wfcr1 (:, :, :) = 0.d0

! start loop over non-reduced k-point set
!      Do ik = 1, nkptnr
       do iq = 1, kqset%nkpt
         ik  = kset%ikp2ik(ikp) ! 1d reduced index -> 1d non-reduced k-point index
         jk  = kqset%kqid(ik,iq) ! k-dependent weight of each q-point???


! determine q-vector
         v (:) = kqset%vkc (:, ik) - kqset%vkc (:, jk)

         Do ig = 1, ngvec
! determine G+q vectors
            vgqc (:, ig) = vgc (:, ig) + v (:) ! Checked: vgc == Gset%vgc

! G+q-vector length and (theta, phi) coordinates
            Call sphcrd (vgqc(:, ig), gqc(ig), tpgqc(:, ig))
! spherical harmonics for G+q-vector
            Call genylm (input%groundstate%lmaxvr, tpgqc(:, ig), ylmgq(:, ig))
         End Do
! structure factors for G+q
         Call gensfacgp (ngvec, vgqc, ngvec, sfacgq)
! find the shortest G+q-vector
         Call findigp0 (ngvec, gqc, igq0)
         sfacgq0 (:) = sfacgq (igq0, :)
! compute the required spherical Bessel functions
         lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
         Call genjlgpr (lmax, gqc, jlgqr)
!         Call genjlgpr (lmax, gc, jlgqr)
         Call genjlgq0r (gqc(igq0), jlgq0r)
! calculate the wavefunctions for occupied states
         call genWF(ik,wf1)
         call genWFinMT(wf1)
         call genWFonMesh(wf1)


         call genWF(jk,wf2)
         call genWFinMT(wf2)
         call genWFonMesh(wf2)
         Do ist2 = 1, nomax !nstsv
!            If (evalfv(ist2,jk) .Lt. efermi) Then
               Do ist3 = 1, nstfv
!                  If (evalsvp(ist2) .Gt. efermi) Then
! calculate the complex overlap density



!-------------------------------------------------------------------
call timesec(ta)
                     call WFprodrs(ist2,wf2,ist3,wf1,prod)
call timesec(tb)
if (ik.eq.jk) then
! write(*,*) 'WFprod',tb-ta
call timesec(ta)
                     Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                    & igq0), sfacgq0, prod%mtrlm(:,:,:,1), prod%ir(:,1), zrho01)
call timesec(tb)

call timesec(ta)
                     prod%ir(:,1)=prod%ir(:,1)-zrho01
                     prod%mtrlm(1,:,:,1)=prod%mtrlm(1,:,:,1)-zrho01/y00

call timesec(tb)
endif

call timesec(ta)
! calculate the Coulomb potential
                     Call zpotcoul (nrcmt, nrcmtmax, nrcmtmax, rcmt, &
                    & igq0, gqc, jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:,:,:,1), &
                    & prod%ir(:,1), pot%mtrlm(:,:,:,1), pot%ir(:,1), zrho02)
call timesec(tb)

if (ik.eq.jk) then
                  Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                  & igq0), sfacgq0, pot%mtrlm(:,:,:,1), pot%ir(:,1), zrho01)

                  pot%ir(:,1)=pot%ir(:,1)-zrho01
                  pot%mtrlm(1,:,:,1)=pot%mtrlm(1,:,:,1)-zrho01/y00
endif
!-------------------------------------------------------------------
                        call genWFonMeshOne(pot)
                       pot%ir=conjg(pot%ir)
                       pot%mtmesh=conjg(pot%mtmesh)
                        call WFprodrs(1,pot,ist2,wf2,prod)

! ------------------------------------------------------------------
! ------------------------------------------------------------------
! ------------------------------------------------------------------
                        zvclir(:,ist3)=zvclir(:,ist3)+prod%ir(:,1)*wkptnr(jk)
                        zvclmt(:,:,:,ist3)=zvclmt(:,:,:,ist3)+prod%mtrlm(:,:,:,1)*wkptnr(jk)
! ------------------------------------------------------------------
! ------------------------------------------------------------------
! ------------------------------------------------------------------



if (.false.) then
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
                     Do ist1 = 1, nstsv

!                        If (evalsvp(ist1) .Lt. efermi) Then
! calculate the complex overlap density
call timesec(ta)
                            call WFprod(ist3,wf2,ist1,wf1,prod)
call timesec(tb)
!write(*,*) 'WFprod',tb-ta
call timesec(ta)
                            Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                           & igq0), sfacgq0, prod%mtrlm(:,:,:,1), prod%ir(:,1), zrho01)
call timesec(tb)
!write(*,*) 'zrhogp',tb-ta
call timesec(ta)
                            prod%ir(:,1)=prod%ir(:,1)-zrho01
                            prod%mtrlm(1,:,:,1)=prod%mtrlm(1,:,:,1)-zrho01/y00
call timesec(tb)
!write(*,*) 'remove average',tb-ta
call timesec(ta)
                           zt1 = zfinp (.True., prod%mtrlm(:,:,:,1), zvclmt, prod%ir(:,1), zvclir)
call timesec(tb)
!write(*,*) 'zfinp',tb-ta

!stop
!-------------------------------------------------------------------
! compute the density coefficient of the smallest G+q-vector
                           zt2 = cfq * wiq2 (iq) * &
                          & (conjg(zrho01)*zrho02)
                           zt2 =0d0
                           vnlvv (ist1, ist2) = vnlvv (ist1, ist2) - &
                          & (wkptnr(ik)*zt1+zt2)
!                        End If
                     End Do
end if



! end loop over ist2
!                  End If
               End Do
! end loop over ist3
!            End If
         End Do
! end loop over non-reduced k-point set
      End Do


Do ist1 = 1, nstsv
      Do ist3 = 1, nstsv
            zt1 = zfinp (.True., wf1%mtrlm(:,:,:,ist1),zvclmt(:,:,:,ist3), wf1%ir(:,ist1), zvclir(:,ist3))
            vnlvv (ist1, ist3) = vnlvv (ist1, ist3) - zt1
      End Do
End Do



!----------------------------------------------!
!     valence-core-valence contribution     !
!----------------------------------------------!
! begin loops over atoms and species
  if (.true.) then
  vnlvv(:,:) = 0.d0
  zvclmt(:, :, :, :) = 0.d0
  wfcr1 = 0.d0

  Do is = 1, nspecies
    nrc = nrcmt(is)
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      Do ist2 = 1, spnst (is)
        If (spcore(ist2, is)) Then
           Do m1 = spk (ist2, is) - 1, spk (ist2, is) - 1
           ! Do m1 = - spk (ist2, is), spk (ist2, is) - 1

  ! pass m-1/2 to wavefcr
              Call wavefcr (input%groundstate%lradstep, is, ia, &
             & ist2, m1, nrcmtmax, wfcr1) !Psi*_{a}; Returns in SC (I think)

            ! shape(wfcr1) = (lmmaxvr, nrcmtmax, 2) = (169, 600, 2)
            ! shape(zfshtvr) = (lmmaxvr, lmmaxvr)
            ! shape(mtmesh) = (ntpll,nrmtmax,natmtot,nstsv) = (169, 600, 1, 52)
            ! shape(zfmt) = (lmmaxvr, nrcmtmax) = (169, 600)
            ! shape(zrhomt) = (lmmaxvr, nrcmtmax, natmtot) = (169, 600, 1)
            ! shape(zvclmt) = (lmmaxvr, nrcmtmax, natmtot, nstsv)
            ! shape(wf%mtrlm) = (lmmaxvr, nrmtmax,natmtot,nstsv)
            ! zbshtvr(ntp,lmmaxvr) x wfcr1(lmmaxvr,nrad) = (169,169)x(169,600)
            ! = (169,600) = (ntp,nrad) = zvcltp

              Do ist3 = 1, nstsv
                 ! If (evalfv(ist3,ik) .Gt. efermi) Then

  ! calculate the complex overlap density
                   Call vnlrhomt (.true., is, wfcr1(:, :, 1), &
                  & wf1%mtmesh(:, :, ias, ist3), zrhomt(:, :, &
                  & ias)) ! Psi*_{a}.Psi_{nk} = rho_{a;nk}; Returns in SH

                  ! Do ist1 = 1,nrc
                  !   Write(*,*) dble(zrhomt(ist1,100,ias)), dble(zrhomt(ist1,100,ias))
                  ! End do
                  ! stop

  ! calculate the Coulomb potential
                    zfmt = zrhomt(:,:,1)
                   !  Call zpotclmt (input%groundstate%ptnucl, &
                   ! & input%groundstate%lmaxvr, nrc, rcmt(:, is), &
                   ! & 0.d0, lmmaxvr, zrhomt(:, :, ias), zfmt) ! Returns SH

                   Call zgemm ('N', 'N', lmmaxvr, nrc, lmmaxvr, &
                   & zone, zbshtvr, lmmaxvr, zfmt, lmmaxvr, &
                   & zzero, zvcltp, lmmaxvr) ! Returns zvcltp in SC

  ! calculate the complex overlap density
                    zfmt = zvcltp
                    Call vnlrhomt (.False., is, wfcr1(:, :, 1), zfmt, &
                   & zrhomt(:, :, ias))

                   zvclmt(:, :, ias, ist3) =  zvclmt(:, :, ias, ist3) + zrhomt(:, :, ias)
  ! end loop over ist3
                 ! End If
              End Do
  ! end loops over ist2 and m1
           End Do
        End If
      End do
    End Do
  End Do

  zt1 = 0.d0
  zvcltp = 0.d0
  Do is = 1, nspecies
    nrc = nrcmt(is)
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      Do ist3 = 1, nstsv

        ! convert muffin-tin potential to spherical coordinates
        ! Call zgemm ('N', 'N', lmmaxvr, nrc, lmmaxvr, &
        ! & zone, zbshtvr, lmmaxvr, zvclmt(:, :, ias, ist3), lmmaxvr, &
        ! & zzero, zvcltp, lmmaxvr)
        zvcltp = zvclmt(:, :, ias, ist3)

        Do ist1 = 1, nstsv
             !  zt1 = zfmtinp (.False., &
             ! & input%groundstate%lmaxvr, nrc, &
             ! & rcmt(:, is), lmmaxvr, &
             ! & wf1%mtrlm(:,:,ias,ist1), zvcltp)
             zt1 = zfmtinp (.False., &
            & input%groundstate%lmaxvr, nrc, &
            & rcmt(:, is), lmmaxvr, &
            & wf1%mtmesh(:,:,ias,ist1), zvcltp)
              vnlvv (ist1, ist3) = vnlvv (ist1, ist3) - zt1
  ! End do of ist1
        End Do
  ! End do of ist3
      End Do
    End do
  End do
  End if


      Deallocate (igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
      Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
!      Deallocate (evalsvp)
      Deallocate (evalsvnr)
!     Deallocate (evecfv, evecsv)
!      Deallocate (apwalm)
      Deallocate (sfacgknr, ylmgq, sfacgq)
      Deallocate (wfmt1, wfmt2, wfir1, wfir2, wfcr1, wfcr2)
      Deallocate (zrhomt, zrhoir, zvclmt, zvclir, zvcltp, zfmt)
      call WFRelease(wf1)
      call WFRelease(wf2)
      call WFRelease(prod)
write(*,*) 'WFRelease done'
      deallocate(gntyyy)
      Return
End Subroutine
!EOC
