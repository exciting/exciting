!BOP
!
!!ROUTINE: calcbarcmb
!
!!INTERFACE:
!
subroutine calcbarcmb(iq)
!
!!DESCRIPTION:
!
!This subroutine calculates the matrix of the bare coulomb potential
!
!!USES:
    use modinput
    use modgw
    use mod_mpi_gw, only: myrank
    
!!INPUT PARAMETERS: 
    implicit none
    integer, intent(in) :: iq ! index of the q-point

!!LOCAL VARIABLES:
    integer :: imix, jmix, igq, jgq
    real(8) :: tstart, tend, t0, t1
    ! for diagonalization subroutine
    real(8) :: vl, vu, abstol
    integer :: il, iu, neval, lwork, info, lrwork, liwork
    complex(8), allocatable :: work(:)
    real(8),    allocatable :: rwork(:)
    integer,    allocatable :: iwork(:), ifail(:), isuppz(:)

    real(8), external :: dlamch
      
!!REVISION HISTORY:
! 
! Created Jan 2014 by DIN
!
!EOP
!BOC      
    call timesec(tstart)
        
!=============================================================================== 
! Setup the bare Coulomb potential matrix in MB representation 
!===============================================================================

    if (allocated(barc)) deallocate(barc)
    allocate(barc(matsiz,matsiz))
    barc(:,:) = 0.d0
    
    select case (trim(input%gw%barecoul%basis))
    
    case('pw')
    
      call calcmpwmix(iq)
      call calcbarcmb_pw(iq)
      
    case('mb')
    
      if (Gamma) then
        !------------------------------------------------
        ! Matrix elements for the singular q=0, L=0 case
        !------------------------------------------------
        !call timesec(t0)
        call barcq0
        !call timesec(t1)
        !write(*,*) 'barcq0', t1-t0 
      end if
        
      !-----------------------------------------------------------
      ! Matrix elements between MT and MT mixed product functions
      !-----------------------------------------------------------
      !call timesec(t0)
      call calcbarcmb_mt_mt(iq)
      !call timesec(t1)
      !write(*,*) 'calcbarcmb_mt_mt', t1-t0 
    
      !-----------------------------------------------------------
      ! Matrix elements between an atomic mixed function and an IPW
      !-----------------------------------------------------------
      !call timesec(t0)
      call calcbarcmb_ipw_mt(iq)
      !call timesec(t1)
      !write(*,*) 'calcbarcmb_ipw_mt', t1-t0 
    
      !-----------------------------------------------------------
      ! Matrix elements between two IPW's
      !-----------------------------------------------------------
      !call timesec(t0)
      call calcbarcmb_ipw_ipw(iq)
      !call timesec(t1)
      !write(*,*) 'calcbarcmb_ipw_ipw', t1-t0
      
    case default
    
      write(*,*) 'ERROR(calcbarcmb): Unknown basis type!'
      stop
      
    end select
    
!===============================================================================
! Diagonalize the bare coulomb matrix
!===============================================================================

    !call timesec(t0)

    if (allocated(vmat)) deallocate(vmat)
    allocate(vmat(matsiz,matsiz))
      
    if (allocated(barcev)) deallocate(barcev)
    allocate(barcev(matsiz))
    
if (.false.) then 
    vmat(1:matsiz,1:matsiz) = barc(1:matsiz,1:matsiz)
    deallocate(barc)  
    lwork = 2*matsiz
    allocate(work(lwork),rwork(3*matsiz))
    call zheev( 'v','u',matsiz,vmat,matsiz, &
    &           barcev,work,lwork,rwork,info)
    call errmsg(info.ne.0,'CALCBARCMB',"Fail to diag. barc by zheev !!!")
    deallocate(work,rwork)
else
    lrwork = -1
    liwork = -1
    lwork = -1
    iu = matsiz
    abstol = 2.d0*dlamch('S')
    allocate(work(1),rwork(1),iwork(1),isuppz(1))
    call zheevr('V', 'A', 'U', matsiz, barc, matsiz, vl, vu, il, iu, &
    &           abstol, neval, barcev, vmat, matsiz, isuppz, work, lwork, rwork, &
    &           lrwork, iwork, liwork, info)
    call errmsg(info.ne.0,'CALCBARCMB',"Fail to diag. barc by zheevr !!!")
    lrwork=int(rwork(1))
    liwork=int(iwork(1))
    lwork=int(work(1))
    ! write(*,*) lrwork,liwork,lwork
    deallocate(work,rwork,iwork,isuppz)
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    allocate(isuppz(2*matsiz))
    call zheevr('V', 'A', 'U', matsiz, barc, matsiz, vl, vu, il, iu, &
    &           abstol, neval, barcev, vmat, matsiz, isuppz, work, lwork, rwork, &
    &           lrwork, iwork, liwork, info)
    call errmsg(info.ne.0,'CALCBARCMB',"Fail to diag. barc by zheevr !!!")
    deallocate(work,rwork,iwork,isuppz)
    deallocate(barc)
end if

    !call timesec(t1)
    !write(*,*) 'barc diagonalization', t1-t0

!----------------------    
! debug info
!----------------------
    if (input%gw%debug) then
      !----------------------    
      ! Memory usage info
      !----------------------
      msize = sizeof(barcev)*b2mb+sizeof(vmat)*b2mb
      write(*,'("calcbarcmb: rank, size(Coulomb potential) (Mb):",i4,f12.2)') myrank, msize
      write(fdebug,*) "### barcev ###"
      do imix = 1, matsiz
        write(fdebug,'(i5,e16.6)') imix, barcev(imix) 
      end do
      write(fdebug,*) "### vmat ###"
      do imix = 1, matsiz, matsiz/10
        do jmix = 1, matsiz, matsiz/10
          write(fdebug,'(2i5,2e16.6)') imix, jmix, vmat(imix,jmix)
        end do
      end do
    endif !debug

    call timesec(tend)
    time_barcmb = time_barcmb+tend-tstart
    
end subroutine ! calcbarcmb
!EOC
