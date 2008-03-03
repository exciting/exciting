
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
!!$    real(8) :: s(3,3)
    integer :: oc1,oc2

! *** old stuff ***
!!$    if (symwings) then
!!$!       s(:,:)=0.5d0*(symdfq0(:,:)+transpose(symdfq0))
!!$       s(:,:)=symdfq0(:,:)
!!$       pout=zzero
!!$       puot=zzero
!!$       do i=1,3
!!$          pout=pout+s(i,i)*pou(i)
!!$          puot=pout+s(i,i)*puo(i)
!!$       end do
!!$    else
!!$       pout=pou(optcomp(1,1))
!!$       puot=puo(optcomp(2,1))
!!$    end if

    oc1=optcomp(1,1)
    oc2=optcomp(2,1)
    if (symwings) then
       write(*,'("Error(dfqoscwg): symmetrization of wings under development")')
       call terminate
    else
       pout1=pou(oc1)
       pout2=pou(oc2)
       puot1=puo(oc1)
       puot2=puo(oc2)
    end if

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
