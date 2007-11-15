
module m_fxc_lrc
  implicit none

contains

  subroutine fxc_lrc(msiz,sw,alpha,fxc)
    !
    ! Static long range xc-kernel.
    ! Calculates either the symmetrized fxc_(G,Gp) = -(alpha/4pi)*delta_(G,Gp),
    ! or fxc_(G,Gp) = -(alpha/4pi)*delta_(G,Gp)*delta(G,0).
    !
    use modmain
    use modxs
    implicit none
    ! arguments
    integer, intent(in) :: msiz
    ! true if all G-components of fxc are to be considered
    logical, intent(in) :: sw
    real(8), intent(in) :: alpha
    complex(8), intent(out) :: fxc(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_lrc'
    real(8) :: t1
    integer :: sh(2),n,ig

    sh=shape(fxc)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout,'(a,2i9,a,i9,a)') 'Error('//trim(thisnam)//'): size of &
            &fxc is to small (required)', sh, '(', msiz, ')'
       call terminate
    end if
    
    fxc(:,:)=(0.d0,0.d0)
    t1=-alpha/fourpi
    if (.not.sw) then
       fxc(1,1)=t1
    else
       do ig=1,msiz
          fxc(ig,ig)=t1
       end do
    end if
    
  end subroutine fxc_lrc

end module m_fxc_lrc
