 
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getlocmt
  implicit none
contains

  ! local orbitals functions
  subroutine getlocmt(iq,ik,apwlm)
    use modmain
    use modmpi
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    complex(8), intent(out) :: apwlm(:,:,:,:)
    ! local variables
    character(*), parameter :: thisnam='getlocmt'
    character(256) :: filextt
    integer :: recl
    real(8) :: vql_(3),vkl_(3),vklt(3),vqlt(3)
    real(8), external :: r3dist
    filextt=filext
    if (iq.eq.0) call genfilextread(task)
    inquire(iolength=recl) vql_,vkl_,apwlm
    call getunit(unit1)
    open(unit1,file='LODLM'//trim(filext),action='read',&
         form='unformatted',status='old',access='direct',recl=recl)
    read(unit1,rec=ik) vql_,vkl_,apwlm
    close(unit1)
    if (iq.eq.0) then
       vklt(:)=vkl0(:,ik)
       vqlt(:)=0.d0
    else
       vklt(:)=vkl(:,ik)
       vqlt(:)=vql(:,iq)
    end if
    ! check consistency
!!$    if ((r3dist(vkl_,vklt).gt.epslat).or.(r3dist(vql_,vqlt).gt.epslat)) then
    if (r3dist(vkl_,vklt).gt.epslat) then
       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
            &APW MT coefficients (current/file): '
       if (procs.gt.1) write(unitout,'(a,i6)') '(parallel) rank', rank
       write(unitout,'(a,i6)') ' q-point index  :', iq
       write(unitout,'(a,i6)') ' k-point index  :', ik
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vql            :', vqlt,',', vql_
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vkl            :', vklt,',', vkl_
       write(unitout,'(a)')    ' file           : LODLM'//trim(filext)
       call terminate
    end if
    filext=filextt
  end subroutine getlocmt

end module m_getlocmt
