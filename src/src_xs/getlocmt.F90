 
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getlocmt
  implicit none
contains

  ! local orbitals functions
  subroutine getlocmt(iq,ik,isti,istf,lolm)
    use modmain
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,isti,istf
    complex(8), intent(out) :: lolm(:,:,:,:)
    ! local variables
    character(*), parameter :: thisnam='getlocmt'
    character(256) :: filextt
    integer :: un,recl,err,nstfv_,nlomax_,lolmax_
    real(8) :: vql_(3),vkl_(3),vklt(3),vqlt(3)
    complex(8), allocatable :: lolmt(:,:,:,:)
    real(8), external :: r3dist
    err=0
    ! check band range
    if ((isti.lt.1).or.(istf.gt.nstfv).or.(istf.le.isti)) then
       write(unitout,*)
       write(unitout,'("Error(getlocmt): inconsistent limits for bands:")')
       write(unitout,'(" band limits  : ",2i6)') isti,istf
       write(unitout,'(" maximum value: ",i6)') nstfv
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if (size(lolm,1).ne.(istf-isti+1)) then
       write(unitout,*)
       write(unitout,'("Error(getlocmt): output array does not match for &
            &bands:")')
       write(unitout,'(" band limits              : ",2i6)') isti,istf
       write(unitout,'(" requested number of bands: ",i6)') istf-isti+1
       write(unitout,'(" array size               : ",i6)') size(lolm,1)
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    ! set file extension
    filextt=filext
    if (iq.eq.0) call genfilextread(task)
    !------------------------!
    !     get parameters     !
    !------------------------!
    inquire(iolength=recl) vql_,vkl_,nstfv_,nlomax_,lolmax_
    call getunit(un)
    open(un,file='LOCMT'//trim(filext),action='read',form='unformatted', &
         status='old',access='direct',recl=recl)
    read(un,rec=1) vql_,vkl_,nstfv_,nlomax_,lolmax
    close(un)
    err=0
    ! check number of bands
    if (nstfv.gt.nstfv_) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): invalid nstfv for k-point ",I8)') &
            thisnam,ik
       write(unitout,'(" q-point    : ",I8)') iq
       write(unitout,'(" current    : ",I8)') nstfv
       write(unitout,'(" FILE       : ",I8)') nstfv_
       write(unitout,'(" filename   : ",a )') 'LOCMT'//trim(filext)
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    ! check number of local orbitals
    if (nlomax.ne.nlomax_) then
       write(unitout,*)
       write(unitout,'("Error(",a,"): invalid nlomax for k-point ",I8)') &
            thisnam,ik
       write(unitout,'(" q-point    : ",I8)') iq
       write(unitout,'(" current    : ",I8)') nlomax
       write(unitout,'(" FILE       : ",I8)') nlomax_
       write(unitout,'(" filename   : ",a )') 'LOCMT'//trim(filext)
       write(unitout,*)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! assign to output array and apply cutoff
    allocate(lolmt(nstfv_,nlomax,-lolmax:lolmax,natmtot))
    ! read data from file
    inquire(iolength=recl) vql_,vkl_,nstfv_,nlomax_,lolmax_,lolmt
    call getunit(un)
    open(un,file='LOCMT'//trim(filext),action='read',form='unformatted', &
         status='old',access='direct',recl=recl)
    read(un,rec=ik) vql_,vkl_,nstfv_,nlomax_,lolmax_,lolmt
    close(un)
    ! check q-point and k-point
    if (iq.eq.0) then
       ! Gamma Q-point
       vklt(:)=vkl0(:,ik)
       vqlt(:)=0.d0
    else
       vklt(:)=vkl(:,ik)
       vqlt(:)=vql(:,iq)
    end if
    if ((r3dist(vkl_,vklt).gt.epslat).or.(r3dist(vql_,vqlt).gt.epslat)) then
       write(unitout,*)
       write(unitout,'(a)') 'Error('//thisnam//'): differring parameters for &
            &LO MT coefficients (current/file): '
       write(unitout,'(a,i6)') ' q-point index  :', iq
       write(unitout,'(a,i6)') ' k-point index  :', ik
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vql            :', vqlt,',', vql_
       write(unitout,'(a,3f12.6,a,3f12.6)') ' vkl            :', vklt,',', vkl_
       write(unitout,'(a)')    ' file           : LOCMT'//trim(filext)
       write(unitout,*)
       call flushifc(unitout)
       call terminate
    end if
    ! retreive data within cutoffs
    lolm(:,:,:,:)=lolmt(isti:istf,:,:,:)
    deallocate(lolmt)
    ! restore file extension
    filext=filextt
  end subroutine getlocmt
end module m_getlocmt
