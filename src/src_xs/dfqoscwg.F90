
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
    complex(8) :: pout,puot
    real(8) :: s(3,3)
    integer :: i

    if (symwings) then
!!$       s(:,:)=0.5d0*(symdfq0(:,:)+transpose(symdfq0))
       s(:,:)=symdfq0(:,:)
       pout=zzero
       puot=zzero
       do i=1,3
          pout=pout+s(i,i)*pou(i)
          puot=pout+s(i,i)*puo(i)
       end do
    else
!!$       oc=optcomp(1,1)
!!$       pout=pou(oc)
!!$       puot=puo(oc)
       pout=pou(optcomp(1,1))
       puot=puo(optcomp(2,1))
    end if

    if (sw.eq.1) then
       ! first wing (G=0,G/=0)
       you(:)=pout*conjg(xou(:))
       yuo(:)=puot*conjg(xuo(:))
    else
       ! second wing (G/=0,G=0)
       you(:)=xou(:)*conjg(pout)
       yuo(:)=xuo(:)*conjg(puot)
    end if

  end subroutine dfqoscwg

end module m_dfqoscwg
