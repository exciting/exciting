
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linkfilesq0
  implicit none
  call system('ln -sf EVECFV_QMT001.OUT  EVECFV_QMT000.OUT')
  call system('ln -sf EVECSV_QMT001.OUT  EVECSV_QMT000.OUT')
  call system('ln -sf EVALSV_QMT001.OUT  EVALSV_QMT000.OUT')
  call system('ln -sf OCCSV_QMT001.OUT   OCCSV_QMT000.OUT')
  call system('ln -sf APWCMT_QMT001.OUT  APWCMT_QMT000.OUT')
  call system('ln -sf LOCMT_QMT001.OUT   LOCMT_QMT000.OUT')
  call system('ln -sf EIGVAL_QMT001.OUT  EIGVAL_QMT000.OUT')
  call system('ln -sf KPOINTS_QMT001.OUT KPOINTS_QMT000.OUT')
  call system('ln -sf EFERMI_QMT001.OUT  EFERMI_QMT000.OUT')
end subroutine linkfilesq0
