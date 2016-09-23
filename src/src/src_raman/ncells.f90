! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted for exciting 2014, Stefan Kontur
!
!
subroutine NCELLS(sn,nc)
!
   use raman_coeff, only : a0,a1,a2,a3,a4,a5,a6
   implicit none
   real(8), intent(in) :: sn
   integer, intent(in) :: nc
!
!
!  determine coefficients for N cells (i.e. the coherence volume)
!  assumes minimum of potential at x=0 (hence a1 can be set to 0)
!
!  nc...number of cells, sn...sqrt( nc )
      write(66,90)
      a0 = 0.0d0 
!     a1 = a1*sn = 0.d0
      a1 = 0.d0
      a2 = a2
      a3 = a3/sn
      a4 = a4/dble(nc)
      a5 = a5/sn/dble(nc)
      a6 = a6/dble(nc)/dble(nc)
      write(66,95) nc,a0,a1,a2,a3,a4,a5,a6
      return
  90  format(//,53('*'),'  NCELLS  ',53('*')/)
  95  format(/' Coeff. for ',I10,' cells:',7f11.5)
!100  format(//,116('-'),//)
end subroutine ncells
!
!
