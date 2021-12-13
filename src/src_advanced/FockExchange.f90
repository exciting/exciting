!
!
!
! Copyright (C) 2021 D. Zavickis, A. Gulans
!
!
Subroutine FockExchange (ikp, vnlvv, vxpsiir,vxpsimt)
      Use modmain
      Use modinput
      Use modgw, only : kqset,Gkqset, kset, nomax, numin, ikvbm, ikcbm, ikvcm, Gset
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (Out) :: vnlvv (nstsv, nstsv)
      Complex (8), Intent (Out) :: vxpsiir (ngrtot, nstsv)
      Complex (8), Intent (Out) :: vxpsimt (lmmaxvr, nrcmtmax, natmtot, nstsv)

! local variables
      Integer :: ngknr, ik, jk, ist1, ist2, ist3
      Integer :: is, ia, ias, ic, m1, m2, lmax, lm
      Integer :: nrc, iq, ig, iv (3), igq0, igk
      Integer :: ilo, loindex
      Integer :: info

      Real (8) :: v (3), cfq, ta,tb
      Complex (8) zrho01, zrho02, zt1, zt2
      Integer :: nr, l, m, io1, lm2, ir, if3
! automatic arrays
      Real (8) :: zn (nspecies)
      Complex (8) sfacgq0 (natmtot)
! allocatable arrays
      Real (8), Allocatable :: vgqc (:, :)
      Real (8), Allocatable :: tpgqc (:, :)
      Real (8), Allocatable :: gqc (:)
      Real (8), Allocatable :: jlgqr (:, :, :)
      Real (8), Allocatable :: jlgq0r (:, :, :)
      Real (8), Allocatable :: evalfv (:,:)
      Complex (8), Allocatable :: ylmgq (:, :)
      Complex (8), Allocatable :: sfacgq (:, :)
      Complex (8), Allocatable :: wfcr1 (:, :, :)
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
      Allocate (vgqc(3, ngvec))
      Allocate (tpgqc(2, ngvec))
      Allocate (gqc(ngvec))
      Allocate (jlgqr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1, ngvec, nspecies))
      Allocate (jlgq0r(0:input%groundstate%lmaxvr, nrcmtmax, nspecies))
      Allocate (ylmgq(lmmaxvr, ngvec))
      Allocate (sfacgq(ngvec, natmtot))
      Allocate (wfcr1(ntpll, nrcmtmax, 2))
      Allocate (zrhomt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvcltp(ntpll, nrcmtmax))
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

      allocate(pot%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
      allocate(pot%ir(ngrtot,1))

! factor for long-range term
      cfq = 0.5d0 * (omega/pi) ** 2
! set the nuclear charges to zero
      zn (:) = 0.d0
      vxpsiir (:, :) = 0.d0
      vxpsimt (:, :, :, :) = 0.d0
      vnlvv (:, :) = 0.d0
! calculate the wavefunctions for all states for the input k-point

! start loop over non-reduced k-point set
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
         Call genjlgq0r (gqc(igq0), jlgq0r)
! calculate the wavefunctions for occupied states
         call genWF(ik,wf1)
         call genWFinMT(wf1)
         call genWFonMesh(wf1)

         call genWF(jk,wf2)
         call genWFinMT(wf2)
         call genWFonMesh(wf2)
         Do ist2 = 1, nomax 
               Do ist3 = 1, nstfv

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
!write(*,*) 'zrhogp',tb-ta
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
                        vxpsiir(:,ist3)=vxpsiir(:,ist3)+prod%ir(:,1)*wkptnr(jk)
                        vxpsimt(:,:,:,ist3)=vxpsimt(:,:,:,ist3)+prod%mtrlm(:,:,:,1)*wkptnr(jk)
! ------------------------------------------------------------------
! ------------------------------------------------------------------
! ------------------------------------------------------------------

! end loop over ist3
               End Do
! end loop over ist2
         End Do
! end loop over non-reduced k-point set
      End Do

!----------------------------------------------!
!     valence-core-valence contribution     !
!----------------------------------------------!
! begin loops over atoms and species

If (.true.) Then
Do is = 1, nspecies
  nrc = nrcmt(is)
  Do ia = 1, natoms (is)
    ias = idxas (ia, is)
    Do ist2 = 1, spnst (is)
      If (spcore(ist2, is)) Then
         Do m1 = - spk (ist2, is), spk (ist2, is) - 1
! pass m-1/2 to wavefcr
            Call wavefcr2 (input%groundstate%lradstep, is, ia, &
           & ist2, m1, nrcmtmax, wfcr1) !Psi*_{a}; Returns in SC (I think)

! Begin loop over occupied and empty states
            Do ist3 = 1, nstsv
               ! If (evalfv(ist3,ik) .Gt. efermi) Then

! calculate the complex overlap density

                Call vnlrhomt2 (.true., is, wfcr1(:, :, 1), &
                & wf1%mtmesh(:, :, ias, ist3), zrhomt(:, :, &
                & ias)) ! Psi*_{a}.Psi_{nk} = rho_{a;nk}; Returns in SH)

! calculate the Coulomb potential

                Call zpotclmt (input%groundstate%ptnucl, &
                & input%groundstate%lmaxvr, nrc, rcmt(:, is), &
                & 0.d0, lmmaxvr, zrhomt(:, :, ias), zfmt) ! Returns SH

                Call zgemm ('N', 'N', ntpll, nrc, lmmaxvr, &
                & zone, zbshthf, ntpll, zfmt, lmmaxvr, &
                & zzero, zvcltp, ntpll) ! Returns zvcltp in SC

! calculate the complex overlap density
                Call vnlrhomt2 (.true., is, wfcr1(:, :, 1), zvcltp, &
                & zrhomt(:, :, ias)) ! Returns in SH

                zvclmt(:,:,ias,ist3)=zvclmt(:,:,ias,ist3)+zrhomt(:, :, ias)
! end loop over ist3
                 ! End If
            End Do
! end loops over ist2 and m1
         End Do
      End If
    End do
  End Do
End Do
End If

Do ist1 = 1, nstsv
      Do ist3 = 1, nstsv
            zt1 = zfinp (.True., wf1%mtrlm(:,:,:,ist1),vxpsimt(:,:,:,ist3), wf1%ir(:,ist1), vxpsiir(:,ist3))
            vnlvv (ist1, ist3) = vnlvv (ist1, ist3) - zt1
      End Do
End Do


      Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
      Deallocate (ylmgq, sfacgq)
      Deallocate (wfcr1)
      Deallocate (zrhomt, zrhoir, zvcltp, zfmt)
      call WFRelease(wf1)
      call WFRelease(wf2)
      call WFRelease(prod)

      Return
End Subroutine
!EOC
