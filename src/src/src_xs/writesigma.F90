!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_writesigma
      Implicit None
Contains
!
!
      Subroutine writesigma (iq, w, sigma, fn)
         Use modxs
         Use m_getunit
         Use m_writevars
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq
         Real (8), Intent (In) :: w (:)
         Complex (8), Intent (In) :: sigma (:)
         Character (*), Intent (In) :: fn
    ! local variables
         Character (*), Parameter :: thisnam = 'writesigma'
         Integer :: n1 (1), n, iw
         If (any(shape(w) .Ne. shape(sigma))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input arr&
           &ays have diffenrent shape'
            Call terminate
         End If
         n1 = shape (w)
         n = n1 (1)
         Call getunit (unit1)
         Open (unit1, File=trim(fn), Action='write')
    ! write data to file
         Write (unit1, '(3g18.10)') (w(iw)*escale, sigma(iw), iw=1, n)
    ! write relevant parameters to file
         Call writevars (unit1, iq, iq)
         Close (unit1)
      End Subroutine writesigma
!
End Module m_writesigma
