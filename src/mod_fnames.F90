


module mod_names
use modinput
! filename for first-variational eigenvectors
character(256) :: filetag_evecfv
data filetag_evecfv / 'EVECFV' /
! filename for second-variational eigenvectors
character(256) :: filetag_evecsv
data filetag_evecsv / 'EVECSV' /
! filename for first-variational eigenvalues
character(256) :: filetag_evalfv
data filetag_evalfv / 'EVALFV' /
! filename for second-variational eigenvalues
character(256) :: filetag_evalsv
data filetag_evalsv / 'EVALSV' /
! filename for second-variational occupation numbers
character(256) :: filetag_occsv
data filetag_occsv / 'OCCSV' /

contains


subroutine revert_names
  implicit none
  filetag_evecfv='EVECFV'
  filetag_evecsv='EVECSV'
  filetag_evalfv='EVALFV'
  filetag_evalsv='EVALSV'
  filetag_occsv='OCCSV'
end subroutine revert_names

end module mod_names
