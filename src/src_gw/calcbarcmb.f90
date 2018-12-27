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
    use mod_coulomb_potential
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

    character(len=256) :: filename
      
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

    if (allocated(vmat)) deallocate(vmat)
    allocate(vmat(matsiz,matsiz))
    vmat(:,:) = barc(:,:)
    deallocate(barc)
    
    if (allocated(barcev)) deallocate(barcev)
    allocate(barcev(matsiz))
    
    lrwork = -1
    liwork = -1
    lwork = -1
    iu = matsiz
    abstol = 2.d0*dlamch('S')
    allocate(work(1),rwork(1),iwork(1))
    call zheevd('V', 'U', matsiz, vmat, matsiz, barcev, work, lwork, rwork, lrwork, iwork, liwork, info)
    call errmsg(info.ne.0, 'CALCBARCMB', "Fail to diag. barc by zheevd !!!")

    lrwork=int(rwork(1))
    liwork=int(iwork(1))
    lwork=int(work(1))
    ! write(*,*) lrwork,liwork,lwork
    deallocate(work,rwork,iwork)

    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    call zheevd('V', 'U', matsiz, vmat, matsiz, barcev, work, lwork, rwork, lrwork, iwork, liwork, info)
    call errmsg(info.ne.0, 'CALCBARCMB', "Fail to diag. barc by zheevd !!!")
    deallocate(work,rwork,iwork)

!----------------------    
! debug info
!----------------------

    if (input%gw%debug) then
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
