!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writepmatasc
      Use modinput, only: input
      Use mod_eigenvalue_occupancy, only: nstsv
      Use modxas, only: ncg
      Use mod_kpoint, only: nkpt, vkl
      Use m_getunit, only: getunit
      Use m_getpmat
      Use m_getpmatxas
      Implicit None
      Complex (8), Allocatable :: pmat (:, :, :)
      Integer :: un, ik, ist1, ist2, oct
  ! initialize global variables
      Call init0
      Call init1
      Call init2
      if (.NOT. input%xs%bse%xas) then
        Allocate (pmat(3, nstsv, nstsv))
      else
        Allocate (pmat(3,ncg,nstsv))
      end if
      Call getunit (un)
      Open (un, File='PMAT_XS_ASC.OUT', Action='write')
      Do ik = 1, nkpt
     ! read momentum matrix elements if required
        if (.NOT. input%xs%bse%xas) then 
          Call getpmat (ik, vkl, 1, nstsv, 1, nstsv, .True., 'PMAT_XS.OU&
          &T', pmat)
          Do ist1 = 1, nstsv
            Do ist2 = 1, nstsv
              Do oct = 1, 3
                Write (un, '(3i8, i4, 4g18.10)') ik, ist1, ist2, oct, &
                 & pmat (oct, ist1, ist2), Abs (pmat(oct, ist1, ist2)), &
                 & Abs (pmat(oct, ist2, ist1)-conjg(pmat(oct, ist1, &
                 & ist2)))
              End Do
            End Do
          End Do
        else
          write(*,*) 'Call to getpmatxas' 
            Call getpmatxas(ik, vkl, 1, ncg,1,nstsv, .True., 'PMAT_XS.OUT',pmat)
          Do ist1 = 1, ncg
            Do ist2 = 1, nstsv
              Do oct = 1, 3
!                Write (un, '(3i8, i4, 4g18.10)') ik, ist1, ist2, oct, &
!                 & pmat (oct, ist1, ist2), Abs (pmat(oct, ist1, ist2)), &
!                 & Abs (pmat(oct, ist2, ist1)-conjg(pmat(oct, ist1, &
!                 & ist2)))
                Write (un, '(3i8, i4, 2g18.10)') ik, ist1, ist2, oct, &
                 & pmat (oct, ist1, ist2)
             End Do
            End Do
          End Do
        end if
      End Do
     Close (un)
      Deallocate (pmat)
End Subroutine writepmatasc
