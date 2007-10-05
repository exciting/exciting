
module m_fxc_lrcd
  implicit none

contains

  subroutine fxc_lrcd(msiz,sw,alpha,beta,w,fxc)
    !
    ! Dynamical long range xc-kernel.
    ! Calculates either the symmetrized 
    ! fxc(G,Gp) = -(alpha+beta*w**2)/4pi*delta_(G,Gp),
    ! or fxc_(G,Gp) = -(alpha+beta*w**2)/4pi*delta_(G,Gp)*delta(G,0).
    !
    use modmain
    use modtddft
    implicit none
    ! arguments
    integer, intent(in) :: msiz
    ! true if all G-components of fxc are to be considered
    logical, intent(in) :: sw
    real(8), intent(in) :: alpha,beta
    complex(8), intent(in) :: w
    complex(8), intent(out) :: fxc(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_lrcd'
    complex(8) :: zt1
    integer :: sh(2),n,ig

    sh=shape(fxc)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout,'(a,2i9,a,i9,a)') 'Error('//trim(thisnam)//'): size of &
            &fxc is to small (required)', sh, '(', msiz, ')'
       call terminate
    end if
    
    fxc(:,:)=(0.d0,0.d0)
    zt1=-(alpha+beta*w**2)/fourpi
    if (.not.sw) then
       fxc(1,1)=zt1
    else
       do ig=1,msiz
          fxc(ig,ig)=zt1
       end do
    end if
    
  end subroutine fxc_lrcd

end module m_fxc_lrcd
