
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dfqoscwg
  implicit none
contains

  subroutine dfqoscwg(sw,pou,puo,xou,xuo,you,yuo)
    use modmain
    use modxs
    implicit none
    ! arguments
    integer, intent(in) :: sw
    complex(8), intent(in) :: pou(3),puo(3),xou(:),xuo(:)
    complex(8), intent(out) :: you(:),yuo(:)
    ! local variables
    complex(8) :: pout1,pout2,puot1,puot2
    integer :: oc1,oc2
    oc1=optcomp(1,1)
    oc2=optcomp(2,1)
    pout1=pou(oc1)
    pout2=pou(oc2)
    puot1=puo(oc1)
    puot2=puo(oc2)
    if (sw.eq.1) then
       ! first wing (G=0,G/=0)
       you(:)=pout1*conjg(xou(:))
       yuo(:)=puot1*conjg(xuo(:))
    else
       ! second wing (G/=0,G=0)
       you(:)=xou(:)*conjg(pout2)
       yuo(:)=xuo(:)*conjg(puot2)
    end if
  end subroutine dfqoscwg

end module m_dfqoscwg
