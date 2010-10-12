
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writephn
  use modmain
  implicit none
! initialise universal variables
  Call init0
  Call init2
  call writephn_list(nphwrt,vqlwrt)
end subroutine


subroutine phononinterpolate
  use modinput
  implicit none
  call interphn(input%phonons%interpolate%ngridq, &
    input%phonons%interpolate%vqloff,.true.,.true.)
end subroutine


subroutine interphn(ngridp,vploff,reducep,tfbz)
  implicit none
  ! arguments
  integer, intent(in) :: ngridp(3)
  real(8), intent(in) :: vploff(3)
  logical, intent(in) :: reducep, tfbz
  ! local variables
  real(8) :: boxl(3,4)
  Integer :: nqpti
  Integer, Allocatable :: ivqi (:, :)
  Integer, Allocatable :: iqmapi (:, :, :)
  Real (8), Allocatable :: vqli (:, :)
  Real (8), Allocatable :: vqci (:, :)
  Real (8), Allocatable :: wqpti (:)
  If (allocated(vqli)) deallocate (vqli)
  Allocate (vqli(3, ngridp(1)*ngridp(2)*ngridp(3)))
  Allocate (ivqi(3, ngridp(1)*ngridp(2)*ngridp(3)))
  Allocate (vqci(3, ngridp(1)*ngridp(2)*ngridp(3)))
  Allocate (wqpti(ngridp(1)*ngridp(2)*ngridp(3)))
  Allocate (iqmapi(0:ngridp(1)-1, 0:ngridp(2)-1, 0:ngridp(3)-1))
  boxl (:, 1) = vploff(:) / dble(ngridp(:))
  boxl (:, 2) = boxl (:, 1)
  boxl (:, 3) = boxl (:, 1)
  boxl (:, 4) = boxl (:, 1)
  boxl (1, 2) = boxl (1, 2) + 1.d0
  boxl (2, 3) = boxl (2, 3) + 1.d0
  boxl (3, 4) = boxl (3, 4) + 1.d0
  ! generate q-point set for interpolation
  Call genppts (reducep, tfbz, ngridp, boxl, &
        & nqpti, iqmapi, ivqi, vqli, vqci, wqpti)
  deallocate(ivqi,vqci,wqpti,iqmapi)
  ! interpolate
  call writephn_list(nqpti,vqli)
end subroutine



Subroutine writephn_list(nppt,vpl)
      Use modmain
      Implicit None
! arguments
      integer, intent(in) :: nppt
      real(8), intent(in) :: vpl(3,nppt)
! local variables
      Integer :: n, iq, i, j, is, ia, ip
! allocatable arrays
      Real (8), Allocatable :: w (:)
      Complex (8), Allocatable :: ev (:, :)
      Complex (8), Allocatable :: dynq (:, :, :)
      Complex (8), Allocatable :: dynp (:, :)
      Complex (8), Allocatable :: dynr (:, :, :)
      n = 3 * natmtot
      Allocate (w(n))
      Allocate (ev(n, n))
      Allocate (dynq(n, n, nqpt))
      Allocate (dynp(n, n))
      Allocate (dynr(n, n, ngridq(1)*ngridq(2)*ngridq(3)))
! read in the dynamical matrices
      Call readdyn (.true.,dynq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! Fourier transform the dynamical matrices to real-space
      Call dynqtor (dynq, dynr)
      Open (50, File='PHONON.OUT', Action='WRITE', Form='FORMATTED')
      Do iq = 1, nppt
         Call dynrtoq (vpl(:, iq), dynr, dynp)
         Call dyndiag (dynp, w, ev)
         Write (50,*)
         Write (50, '(I6, 3G18.10, " : q-point, vpl")') iq, vpl &
        & (:, iq)
         Do j = 1, n
            Write (50,*)
            Write (50, '(I6, G18.10, " : mode, frequency")') j, w (j)
            i = 0
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  Do ip = 1, 3
                     i = i + 1
                     If (i .Eq. 1) Then
                        Write (50, '(3I4, 2G18.10, " : species, atom, p&
                       &olarisation, eigenvector")') is, ia, ip, ev (i, &
                       & j)
                     Else
                        Write (50, '(3I4, 2G18.10)') is, ia, ip, ev (i, &
                       & j)
                     End If
                  End Do
               End Do
            End Do
         End Do
         Write (50,*)
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writephn): phonon frequencies and eigenvectors &
     &written to PHONON.OUT")')
      Write (*, '(" for all q-vectors in the phwrite list")')
      Write (*,*)
      Deallocate (w, ev, dynq, dynp, dynr)
      Return
End Subroutine
