
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writesigma
  implicit none
contains

  subroutine writesigma(iq,w,sigma,fn)
    use modxs
    use m_getunit
    use m_tdwriteh
    implicit none
    ! arguments
    integer, intent(in) :: iq
    real(8), intent(in) :: w(:)
    complex(8), intent(in) :: sigma(:)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam='writesigma'
    integer :: n1(1),n,iw
    if (any(shape(w).ne.shape(sigma))) then
       write(unitout,'(a)') 'Error('//thisnam//'): input arrays have &
            &diffenrent shape'
       call terminate
    end if
    n1=shape(w); n=n1(1)
    call getunit(unit1)
    open(unit1,file=trim(fn),action='write')
    ! write relevant parameters as header to file
    call tdwriteh(unit1,iq)
    ! write data to file
    write(unit1,'(3g18.10)') (w(iw)*escale,sigma(iw),iw=1,n)
    close(unit1)
  end subroutine writesigma

end module m_writesigma
