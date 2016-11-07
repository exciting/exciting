!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_writekerr
      Implicit None
Contains
!
!
      Subroutine writekerr (iq, w, kerr, fn)
         use modmpi
         Use modxs
         Use m_getunit
         Use m_writevars
         Use mod_constants, Only: pi
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq
         Real (8), Intent (In) :: w (:)
         Complex (8), Intent (In) :: kerr (:)
         Character (*), Intent (In) :: fn
    ! local variables
         Character (*), Parameter :: thisnam = 'writekerr'
         Integer :: n1 (1), n, iw
         If (any(shape(w) .Ne. shape(kerr))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input arr&
           &ays have diffenrent shape'
            Call terminate
         End If
         n1 = shape (w)
         n = n1 (1)
         Call getunit (unit1)
         Open (unit1, File=trim(fn), Action='write')
    ! write data to file
         !Write (unit1, '(3g18.10)') (w(iw)*escale, kerr(iw)*180.d0/pi, iw=1, n)
         Write (unit1, '(3g18.10)') (w(iw)*escale, kerr(iw)*180.d0/pi, iw=2, n)
    ! write relevant parameters to file
         Call writevars (unit1, iq, iq)
         Close (unit1)
      End Subroutine writekerr
!
End Module m_writekerr
