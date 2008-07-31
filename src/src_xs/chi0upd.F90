
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_chi0upd
  implicit none
contains

!BOP
! !ROUTINE: chi0upd
! !INTERFACE:
  subroutine chi0upd(n,wou,wuo,hou,huo,chi0)
! !INPUT/OUTPUT PARAMETERS:
!   n       : size of Fourier basis (in,integer)
!   wou     : inverse energy denominator (occ to unocc) (in,complex)
!   wuo     : inverse energy denominator (unocc to occ) (in,complex)
!   hou     : oscillators (occ to unocc) (in,complex(n,n))
!   huo     : oscillators (unocc to occ) (in,complex(n,n))
!   chi0    : updated Kohn Sham response function (out,complex(n,n))
! !DESCRIPTION:
!   Updates the Kohn-Sham response function by adding a transition which
!   corresponds to a pair of states $o$ and $u$ and a ${\bf k}$-point:
!   $$ \chi^0_{\bf{GG'}}(\omega) \rightarrow
!      \chi^0_{\bf{GG'}}(\omega) + 
!      h^{ou}_{\bf{GG'}}({\bf k}) w_{\rm ou}(\omega) +
!      h^{uo}_{\bf{GG'}}({\bf k}) w_{\rm uo}(\omega). $$
!   Here, $\chi^0_{\bf{GG'}}$ is the Fourier component of the Kohn-Sham 
!   response function, $h^{ou}_{\bf{GG'}},h^{uo}_{\bf{GG'}}$
!   are the oscillators (products of plane wave and/or momentum matrix elements)
!   for an $o \rightarrow u$ and $u \rightarrow o$ transition whereas
!   $w_{\rm ou},w_{\rm uo}$ are the inverse energy denominators carrying the
!   $\omega$-dependence. The ${\bf q}$-point is assumed to be fixed.
!
! !REVISION HISTORY:
!   Created January 2005 (Sagmeister)
!   Added documentation, July 2008 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: wou,wuo,hou(n,n),huo(n,n)
    complex(8), intent(inout) :: chi0(n,n)
    ! local variables
    integer :: m
    do m=1,n
       call zaxpy(n,wou,hou(1,m),1,chi0(1,m),1)
       call zaxpy(n,wuo,huo(1,m),1,chi0(1,m),1)
    end do
  end subroutine chi0upd

end module m_chi0upd
!EOC
