
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
Module mod_names
      Use modinput
! filename for first-variational eigenvectors
      Character (256) :: filetag_evecfv
      Data filetag_evecfv / 'EVECFV' /
! filename for second-variational eigenvectors
      Character (256) :: filetag_evecsv
      Data filetag_evecsv / 'EVECSV' /
! filename for first-variational eigenvalues
      Character (256) :: filetag_evalfv
      Data filetag_evalfv / 'EVALFV' /
! filename for second-variational eigenvalues
      Character (256) :: filetag_evalsv
      Data filetag_evalsv / 'EVALSV' /
! filename for second-variational occupation numbers
      Character (256) :: filetag_occsv
      Data filetag_occsv / 'OCCSV' /
!
Contains
!
!
      Subroutine revert_names
         Implicit None
         filetag_evecfv = 'EVECFV'
         filetag_evecsv = 'EVECSV'
         filetag_evalfv = 'EVALFV'
         filetag_evalsv = 'EVALSV'
         filetag_occsv = 'OCCSV'
      End Subroutine revert_names
!
End Module mod_names
!
