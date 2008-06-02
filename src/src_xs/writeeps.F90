
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writeeps
  implicit none
contains

  subroutine writeeps(iq,w,eps,fn)
    use modmain
    use modxs
    use m_getunit
    use m_writevars
    implicit none
    ! arguments
    integer, intent(in) :: iq
    real(8), intent(in) :: w(:)
    complex(8), intent(in) :: eps(:)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam='writeeps'
    integer :: n1(1),n,iw
    real(8), allocatable :: imeps(:),kkeps(:)
    if (any(shape(w).ne.shape(eps))) then
       write(unitout,'(a)') 'Error('//thisnam//'): input arrays have &
            &diffenrent shape'
       call terminate
    end if
    n1=shape(w); n=n1(1)
    allocate(imeps(n),kkeps(n))
    ! Kramers-Kronig transform imaginary part
    imeps(:)=aimag(eps(:))
    call kramkron(optcomp(1,1),optcomp(2,1),1.d-8,n,w,imeps,kkeps)
    call getunit(unit1)
    open(unit1,file=trim(fn),action='write')
    write(unit1,'(4g18.10)') (w(iw)*escale,eps(iw),kkeps(iw),iw=1,n)
    ! write relevant parameters to file
    call writevars(unit1,iq)
    close(unit1)
    deallocate(imeps,kkeps)
  end subroutine writeeps

end module m_writeeps
