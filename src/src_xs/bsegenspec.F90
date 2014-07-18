!
!
!
! Copyright (C) 2014 S. Kontur and C. Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: bsegenspec
! !INTERFACE:
!
!
Subroutine bsegenspec ()
! !USES:
   Use modmain
   Use modxs
   use modinput
   Use m_genfilname
   Use m_getunit
   Use m_genwgrid
   Use m_genloss
   Use m_gensigma
   Use m_gensumrls
   Use m_writeeps
   Use m_writeloss
   Use m_writesigma
   Use m_writesumrls
! !DESCRIPTION:
!   Collects BSE results from files EXCITON*.OUT and computes
!   spectra anew.
!
! !REVISION HISTORY:
!   Created April 2014, S. Kontur
!EOP
!BOC
   Implicit None
  ! arguments
  ! local variables
   Integer, Parameter :: noptcmp = 3
   integer :: hamsiz, nexc, ivck
   integer :: oct1, oct2, octl, octu, optcompt(3)
   integer :: iw, idum
   real (8) :: sumrls(3), rdum, real_p, imag_p
   character (256) :: fnexc, dotext
   integer :: unexc
   Real (8), Allocatable :: w (:), loss(:, :, :), bse_en(:)
   Complex (8), Allocatable :: oszs (:, :), spectr (:), sigma(:), buf(:,:,:)
  ! external functions
   Integer, External :: l2int
 
   call init0
   call init1
   call init2
  ! size of BSE-Hamiltonian
   hamsiz = (input%xs%bse%nstlbse(2) - input%xs%bse%nstlbse(1) + 1) * &
     &      (input%xs%bse%nstlbse(4) - input%xs%bse%nstlbse(3) + 1) * nkptnr
   nexc =  hamsiz

   Allocate (oszs(nexc, 3), bse_en(nexc))
   Allocate (w(input%xs%energywindow%points), spectr(input%xs%energywindow%points))
   Allocate (buf(3,3,input%xs%energywindow%points))
   Allocate (loss(3, 3, input%xs%energywindow%points), sigma(input%xs%energywindow%points))
   Call genwgrid (input%xs%energywindow%points, input%xs%energywindow%intv, &
     &     input%xs%tddft%acont, 0.d0, w_real=w)
   buf = zzero


    do oct1 = 1, noptcmp
      Call genfilname (basename='EXCITON', tq0=.True., oc1=oct1, &
  &    oc2=oct1, bsetype=input%xs%bse%bsetype, &
  &    scrtype=input%xs%screening%screentype, nar= .Not. &
  &    input%xs%bse%aresbse, filnam=fnexc)
  ! read oscillator strengths
      Call getunit (unexc)
      Open (unexc, File=fnexc, Form='formatted', Action='read', Status='old')
      Do ivck = 1, hamsiz
         Read (unexc, '(i8, 5g18.10)') idum, bse_en(ivck), rdum, rdum, real_p, imag_p
         oszs(ivck, oct1) = cmplx(real_p, imag_p, 8)
         bse_en(ivck) = bse_en(ivck) / escale
         write(*,*) 'exciton ',idum,' : ',bse_en(ivck)
      End Do
      Close (unexc)
    enddo

      Do oct1 = 1, noptcmp
       If (input%xs%dfoffdiag) Then
             octl = 1
             octu = noptcmp
       Else
             octl = oct1
             octu = oct1
       End If
        Do oct2 = octl, octu
         optcompt = (/ oct1, oct2, 0 /)
         spectr (:) = zzero
         Do iw = 1, input%xs%energywindow%points
            Do ivck = 1, nexc
           ! Lorentzian lineshape
               spectr (iw) = spectr (iw) + &
              & oszs(ivck, oct1) * conjg(oszs(ivck, oct2)) * &
              & (1.d0 / (w(iw) - bse_en(ivck) + zi*input%xs%broad) )
            If (input%xs%bse%aresbse) spectr (iw) = spectr (iw) + &
              & oszs(ivck, oct1) * conjg(oszs(ivck, oct2)) * &
              & (1.d0 / (-w(iw) - bse_en(ivck) + zi*input%xs%broad) )
            End Do
         End Do
         spectr (:) = l2int (oct1 .Eq. oct2) * 1.d0 - spectr (:) * 8.d0 * &
        &             pi / omega / nkptnr
         buf(oct1,oct2,:) = spectr(:)
     ! end loops over optical components
        enddo
      End Do
!
      call genloss (buf, loss, noptcmp)
!
   Do oct1 = 1, noptcmp
    If (input%xs%dfoffdiag) Then
          octl = 1
          octu = noptcmp
    Else
          octl = oct1
          octu = oct1
    End If
     Do oct2 = octl, octu
      optcompt = (/ oct1, oct2, 0 /)
      Call genfilname (basename='EPSILON', tq0=.True., oc1=oct1, &
     & oc2=oct2, bsetype=input%xs%bse%bsetype, &
     & scrtype=input%xs%screening%screentype, nar= .Not. &
     & input%xs%bse%aresbse, filnam=fneps)
      Call genfilname (basename='LOSS', tq0=.True., oc1=oct1, &
     & oc2=oct2, bsetype=input%xs%bse%bsetype, &
     & scrtype=input%xs%screening%screentype, nar= .Not. &
     & input%xs%bse%aresbse, filnam=fnloss)
      Call genfilname (basename='SIGMA', tq0=.True., oc1=oct1, &
     & oc2=oct2, bsetype=input%xs%bse%bsetype, &
     & scrtype=input%xs%screening%screentype, nar= .Not. &
     & input%xs%bse%aresbse, filnam=fnsigma)
      Call genfilname (basename='SUMRULES', tq0=.True., oc1=oct1, &
     & oc2=oct2, bsetype=input%xs%bse%bsetype, &
     & scrtype=input%xs%screening%screentype, nar= .Not. &
     & input%xs%bse%aresbse, filnam=fnsumrules)
  ! symmetrize the macroscopic dielectric function tensor
      Call symt2app (oct1, oct2, input%xs%energywindow%points, symt2, buf, spectr)
  ! generate optical functions
      Call gensigma (w, spectr, optcompt, sigma)
      Call gensumrls (w, spectr, sumrls)
  ! write optical functions to file
      Call writeeps (1, oct1, oct2, w, spectr, trim(fneps))
      Call writeloss (1, w, loss(oct1, oct2, :), trim(fnloss))
      Call writesigma (1, w, sigma, trim(fnsigma))
      Call writesumrls (1, sumrls, trim(fnsumrules))
     enddo
   end do

   deallocate (oszs, w, spectr, loss, sigma, buf)

End Subroutine bsegenspec

!EOC
