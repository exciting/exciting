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
    use mod_coulomb_potential
    use modinput
    use modgw
    use modmain
    use modmpi, only: rank
    use precision, only: i32, dp, str_256
#include "offload.fpp"
!!INPUT PARAMETERS:
    implicit none
    integer(i32), intent(in) :: iq ! index of the q-point

!!LOCAL VARIABLES:
    integer(i32) :: imix, jmix, igq, jgq
    real(dp) :: tstart, tend, t0, t1
    character(str_256) :: filename
    complex(dp), allocatable :: barc_lr(:,:)
    ! for diagonalization subroutine
    real(dp) :: vl, vu, abstol
    integer :: il, iu, neval, lwork, info, lrwork, liwork
    complex(dp), allocatable :: work(:)
    real(dp),    allocatable :: rwork(:)
    integer(i32),    allocatable :: iwork(:), ifail(:), isuppz(:)
    real(dp), external :: dlamch

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

    if (allocated(barc)) then
        OMP_OFFLOAD target exit data map(delete: barc)
        deallocate(barc)
    end if
    allocate(barc(matsiz,matsiz), source = zzero)

    select case (trim(input%gw%barecoul%basis))

    case('pw')

      call calcmpwmix(iq)
      call calcbarcmb_pw(iq)

    case('mb')

      if (Gamma) then
        !------------------------------------------------
        ! Matrix elements for the singular q=0, L=0 case
        !------------------------------------------------
        call barcq0
      end if

      !-----------------------------------------------------------
      ! Matrix elements between MT and MT mixed product functions
      !-----------------------------------------------------------
      call calcbarcmb_mt_mt(iq)

      !-----------------------------------------------------------
      ! Matrix elements between an atomic mixed function and an IPW
      !-----------------------------------------------------------
      call calcbarcmb_ipw_mt(iq)

      !-----------------------------------------------------------
      ! Matrix elements between two IPW's
      !-----------------------------------------------------------
      call calcbarcmb_ipw_ipw(iq)

      !----------------------------------------
      ! Only in case of HSE self-consitent run
      !----------------------------------------
      if ((task == 7) .and. (xctype(1)==408)) then
        ! Compute the short-range Coulomb potential
        allocate(barc_lr(matsiz,matsiz))
        call calcbarcmb_lr(iq, barc_lr)
        barc(:,:) = barc(:,:) - barc_lr(:,:)
        deallocate(barc_lr)
      endif

    case default

      write(*,*) 'ERROR(calcbarcmb): Unknown basis type!'
      stop

    end select


!===============================================================================
! Diagonalize the bare coulomb matrix
!===============================================================================

    if (allocated(vmat)) deallocate(vmat)
    call move_alloc(barc, vmat)

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
      write(fdebug,'("calcbarcmb: rank, size(Coulomb potential) (Mb):",i4,f12.2)') rank, msize
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

end subroutine
!EOC
