
module m_writeexciton
  implicit none
contains

  subroutine writeexciton(iq,oct,w,mdf,fn)
    use modmain
    use modxs
    use m_getunit
    use m_tdwriteh
    implicit none
    ! arguments
    integer, intent(in) :: iq,oct
    real(8), intent(in) :: w(:)
    complex(8), intent(in) :: mdf(:)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam = 'writeexcit'
    character(256) ::str,filnam
    real(8), allocatable :: eps(:)
    complex(8), allocatable :: epstet(:)
    real(8) :: wp,we,ep, f,fp,g,gp
    integer :: j,iw,ne,recl

    allocate(eps(nwdos),epstet(nwdos))
    call getunit(unit1)
    open(unit1,file=trim(fn),action='write')
    ! write parameters as header to file
    call tdwriteh(unit1,iq)
    ! write data to file

    do j=1,nexcit(oct)
       wp=excite(j,oct)
       write(unit1,'(3g18.10)') j, wp*escale, excito(j,oct)
    end do

    ! close files
    close(unit1)

    eps(:)=0.d0
    do ne=1,nexcit(oct)
       we=excite(ne,oct)
       do iw=1,nwdf
          ! single Lorentz peak at zero frequency
          eps(iw)=eps(iw)&
               +excito(ne,oct)*1.0d0/pi*brdtd/(brdtd**2+(w(iw)-we)**2) &
               -excito(ne,oct)*1.0d0/pi*brdtd/(brdtd**2+(-w(iw)-we)**2)
          ! convolute with two Lorentz peaks (+/-) at zero frequency
!!$          eps(iw)=eps(iw)+(1.d0/(2*atan(we/brdtd)))*( &
!!$                      brdtd/((w(iw)-we)**2+brdtd**2) - &
!!$                      brdtd/((-w(iw)-we)**2+brdtd**2) )
       end do ! iw
    end do

    ! add continuum part of macr. dielectric function to excitonic part
    eps=eps+aimag(mdf)

    do iw=1,nwdos
       ep=mdfrpa(iw,oct,1)
       f=fxc0(iw,oct)
       fp=fxc0d(iw,oct)
       g=mdfrpa(iw,oct,1)
       gp=mdfrpad(iw,oct)
       write(5000+oct,'(7g18.10)') w(iw)*escale,eps(iw),ep,&
            fxc0(iw,oct),1.d0-1.d0/fxc0(iw,oct),-1.d0/(ep-1.d0),&
            1.d0+fxc0(iw,oct)*(ep-1.d0)
       write(4000+oct,'(7g18.10)') w(iw)*escale,&
            -pi/(f*abs(fp*(g-1.d0)+f*gp)),&
            f,fp,g,gp
    end do
    deallocate(eps,epstet)

  end subroutine writeexciton

end module m_writeexciton
