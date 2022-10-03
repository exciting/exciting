!> Generate the overlap matrix of the product functions in the interstitial region and diagonalize it. 
!>
!> In the interstitial region, the product basis is constructed from planewaves. 
!> These interstitial planewaves (IPW’s) are not orthogonal. To obtain a set of 
!> orthogonal wavefunctions, we diagonalize the overlap matrix by solving the MT
!> eigenvalue equation.
!>
!> The output is the matrix $S_{\vec{G}i}$.
subroutine diagsgi(iq)
  
    use modinput
    use modmain, only : cfunig
    use modgw,   only : Gset, Gqset, sgi, fdebug, time_diagsgi
    implicit none
    integer(4), intent(in) :: iq      

    integer(4) :: ig, ngq
    integer(4) :: igq                  !! Counter: Runs over igq's
    integer(4) :: jgq                  !! Counter: Runs over igq's
    integer(4), dimension(3) :: iig    !! integer coordinates of G-G'
    real(8)   :: tstart, tend
    
    complex(8) :: cfact
    real(8), allocatable :: epsipw(:)
    complex(8), allocatable :: zmat(:,:)
    
    real(8) :: vl, vu, abstol
    integer :: il, iu, neval, lwork, info, lrwork, liwork
    complex(8), allocatable :: work(:)
    real(8),    allocatable :: rwork(:)
    integer,    allocatable :: iwork(:), ifail(:), isuppz(:)
    real(8), external :: dlamch

    external zheev  

    call timesec(tstart)
    
    ngq = Gqset%ngk(1,iq)
    if(allocated(sgi)) deallocate(sgi)
    allocate(sgi(ngq,ngq))
    
    !-----------------------------------------------------------
    ! Calculate the overlap matrix between product plane waves
    !-----------------------------------------------------------
    sgi = 0.0d0
    do igq = 1, ngq
      ig = Gset%ivgig(0,0,0)
      sgi(igq,igq) = conjg(cfunig(ig))
      ! Non diagonal elements:
      do jgq = igq+1, ngq
        iig(:) = Gset%ivg(:,Gqset%igkig(igq,1,iq))- &
        &        Gset%ivg(:,Gqset%igkig(jgq,1,iq))
        ig = Gset%ivgig(iig(1),iig(2),iig(3))
        ! if (ig>??) stop 'G basis is too small for the selected G_{max}^{MB}'
        sgi(igq,jgq) = conjg(cfunig(ig))
        ! NOTE(ALEX) The lower triangle should not be needed by LAPACK CALL - confirm
        !sgi(jgq,igq) = conjg(sgi(igq,jgq))
      end do ! jgq
    end do ! igq
      
    if (input%gw%debug) then
      write(fdebug,*) "diagsgi: iq, ngq=", iq, ngq
      write(fdebug,*) "### sgi-0 ###"
      do jgq = 1, ngq, ngq/10
        do igq = 1, ngq, ngq/10
          write(fdebug,'(2i5,4f12.6)') igq, jgq, sgi(igq,jgq)
        end do
      end do
    end if 

    !-----------------------------------------------------------
    ! Diagonalize sgi
    !-----------------------------------------------------------
    allocate(epsipw(ngq))
    
    lwork = 2 * ngq - 1
    allocate(work(lwork),rwork(3 * ngq - 2))
    call zheev('V','U', ngq, sgi, ngq, epsipw, work, lwork, rwork, info)
    if (info.ne.0) stop "diagsgi: Fail in calling zheev"

    if (input%gw%debug) then
      write(fdebug,*) "### sgi-1 ###"
      do jgq = 1, ngq, ngq/10
        do igq = 1, ngq, ngq/10
          write(fdebug,'(2i5,4f12.6)') igq, jgq, sgi(igq,jgq)
        end do
      end do
      write(fdebug,*)
      write(fdebug,*) "### epsipw ###"
      do igq = 1, ngq
        write(fdebug,'(i5,2f12.6)') igq, epsipw(igq)
      end do 
    end if
    
    !-----------------------------------------------------------
    ! Normalize sgi
    !-----------------------------------------------------------
    do igq = 1, ngq
      cfact = cmplx(1.0d0/sqrt(dabs(epsipw(igq))),0.0d0,8)
      sgi(:,igq) = cfact*sgi(:,igq)
    enddo

    call timesec(tend)
    time_diagsgi = time_diagsgi+tend-tstart
      
end subroutine
