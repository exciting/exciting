 
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getapwcmt
  implicit none
contains
  
  ! APW functions
  subroutine getapwcmt(iq,ik,isti,istf,lmax,apwlm)
    use modmain
    use modmpi
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,isti,istf,lmax
    complex(8), intent(out) :: apwlm(:,:,:,:)
    ! local variables
    character(*), parameter :: thisnam='getapwcmt'
    character(256) :: filextt
    integer :: un,recl,err,nstfv_,apwordmax_,lmaxapw_
    real(8) :: vql_(3),vkl_(3),vklt(3),vqlt(3)
    complex(8), allocatable :: apwlmt(:,:,:,:)
    real(8), external :: r3dist
    err=0
    ! check band range
    if ((isti.lt.1).or.(istf.gt.nstfv).or.(istf.le.isti)) then
       write(*,*)
       write(*,'("Error(getapwcmt): inconsistent limits for bands:")')
       write(*,'(" band limits  : ",2i6)') isti,istf
       write(*,'(" maximum value: ",i6)') nstfv
       write(*,*)
       err=err+1
    end if
    if (size(apwlm,1).ne.(istf-isti+1)) then
       write(*,*)
       write(*,'("Error(getapwcmt): output array does not match for bands:")')
       write(*,'(" band limits              : ",2i6)') isti,istf
       write(*,'(" requested number of bands: ",i6)') istf-isti+1
       write(*,'(" array size               : ",i6)') size(apwlm,1)
       write(*,*)
       err=err+1
    end if
    ! check lmax value
    if ((lmax.gt.lmaxapw).or.(lmax.lt.0)) then
       write(*,*)
       write(*,'(a,i8)') 'Error('//thisnam//'): lmax > lmaxapw or < 0:', &
            lmax
       write(*,*)
       err=err+1
    end if
    if (err.gt.0) call terminate
    ! set file extension
    filextt=filext
    if (iq.eq.0) call genfilextread(task)
    !------------------------!
    !     get parameters     !
    !------------------------!
    inquire(iolength=recl) vql_,vkl_,nstfv_,apwordmax_,lmaxapw_
    call getunit(un)
    open(un,file='APWCMT'//trim(filext),action='read',&
         form='unformatted',status='old',access='direct',recl=recl)
    read(un,rec=1) vql_,vkl_,apwlmt
    close(un)
    err=0
    ! check number of bands
    if (nstfv.gt.nstfv_) then
       write(*,*)
       write(*,'("Error(",a,"): invalid nstfv for k-point ",I8)') thisnam,ik
       write(*,'(" q-point    : ",I8)') iq
       write(*,'(" current    : ",I8)') nstfv
       write(*,'(" FILE       : ",I8)') nstfv_
       write(*,'(" filename   : ",a      )') 'APWCMT'//trim(filext)
       write(*,*)
       err=err+1
    end if
    ! check APW matching order
    if (apwordmax.ne.apwordmax_) then
       write(*,*)
       write(*,'("Error(",a,"): invalid apwordmax for k-point ",I8)') thisnam,ik
       write(*,'(" q-point    : ",I8)') iq
       write(*,'(" current    : ",I8)') apwordmax
       write(*,'(" FILE       : ",I8)') apwordmax_
       write(*,'(" filename   : ",a      )') 'APWCMT'//trim(filext)
       write(*,*)
       err=err+1
    end if
    ! check lmax
    if (lmaxapw.gt.lmaxapw_) then
       write(*,*)
       write(*,'("Error(",a,"): invalid lmaxapw for k-point ",I8)') thisnam,ik
       write(*,'(" q-point    : ",I8)') iq
       write(*,'(" current    : ",I8)') lmaxapw
       write(*,'(" FILE       : ",I8)') lmaxapw_
       write(*,'(" filename   : ",a      )') 'APWCMT'//trim(filext)
       write(*,*)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! assign to output array and apply cutoff
    allocate(apwlmt(nstfv_,apwordmax,(lmaxapw_+1)**2,natmtot))
    ! read data from file
    inquire(iolength=recl) vql_,vkl_,nstfv_,apwordmax_,lmaxapw_,apwlmt
    call getunit(un)
    open(un,file='APWCMT'//trim(filext),action='read',form='unformatted', &
         status='old',access='direct',recl=recl)
    read(un,rec=ik) vql_,vkl_,nstfv_,apwordmax_,lmaxapw_,apwlmt
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
       write(*,'(a)') 'Error('//thisnam//'): differring parameters for &
            &APW MT coefficients (current/file): '
       if (procs.gt.1) write(*,'(a,i6)') '(parallel) rank', rank
       write(*,'(a,i6)') ' q-point index  :', iq
       write(*,'(a,i6)') ' k-point index  :', ik
       write(*,'(a,3f12.6,a,3f12.6)') ' vql            :', vqlt,',', vql_
       write(*,'(a,3f12.6,a,3f12.6)') ' vkl            :', vklt,',', vkl_
       write(*,'(a)')    ' file           : APWCMT'//trim(filext)
       call terminate
    end if
    ! retreive data within cutoffs
    apwlm(:,:,:,:)=apwlmt(isti:istf,:,1:(lmax+1)**2,:)
    deallocate(apwlmt)
    ! restore file extension
    filext=filextt
  end subroutine getapwcmt
end module m_getapwcmt
