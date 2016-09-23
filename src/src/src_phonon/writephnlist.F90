
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine writephnlist(nppt,vpl,twrev,fname)
      Use modmain
      Implicit None
! arguments
      integer, intent(in) :: nppt
      real(8), intent(in) :: vpl(3,nppt)
      logical, intent(in) :: twrev
      character(*), intent(in) :: fname
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
      Open (50, File=trim(fname), Action='WRITE', Form='FORMATTED')
      Do iq = 1, nppt
         Call dynrtoq (vpl(:, iq), dynr, dynp)
         Call dyndiag (dynp, w, ev)
         Write (50,*)
         Write (50, '(I6, 3G18.10, " : q-point, vpl")') iq, vpl &
        & (:, iq)
         Do j = 1, n
            if (twrev) Write (50,*)
            Write (50, '(I6, G18.10, " : mode, frequency")') j, w (j)
            if (twrev) then
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
            end if
         End Do
         Write (50,*)
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writephn): phonon frequencies (and eigenvectors) &
     &written to ",a)') trim(fname)
      Write (*, '(" for all q-vectors in the phwrite list")')
      Write (*,*)
      Deallocate (w, ev, dynq, dynp, dynr)
      Return
End Subroutine
