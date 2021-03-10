! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine write_wfplot
  use modinput, only: input
  use modxs, only: unitout, iqmtgamma
  use wfplot_hdf5, only: write_wfplot_hdf5, write_box_hdf5
  use m_genfilname
  Implicit none
  integer :: i, k

  call genfilname(iqmt=iqmtgamma, setfilext=.true.)

#ifdef _HDF5_
    call write_box_hdf5(input%xs%scrwfplot%plot3d)
#endif
    do i=input%xs%scrwfplot%bandrange(1), input%xs%scrwfplot%bandrange(2) 
#ifdef _HDF5_
      do k=input%xs%scrwfplot%kptrange(1), input%xs%scrwfplot%kptrange(2)
        call write_wfplot_hdf5(k,i, plot3d_=input%xs%scrwfplot%plot3d)
      end do
#endif 
    end do

  write(unitout, '("Info(scr_wfplot_hdf5):&
    & Real-space plot of screening wavefunctions finished")')

end subroutine
