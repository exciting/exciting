!BOP
!
!!ROUTINE: calcinveps
!
!!INTERFACE:
!
subroutine calcinveps(iomstart,iomend)
!
!!DESCRIPTION:
!
! This subroutine calculates the inverse of the dielectric function
!
!!USES:
    use modmain
    use modgw
    use mod_mpi_gw, only : myrank

!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iomstart, iomend

!!LOCAL VARIABLES:
    integer(4) :: i, j, iom
    integer(4) :: im, jm
    integer(4) :: info, lwork
   
    real(8)    :: tstart, tend
    real(8)    :: q0eps(3), modq0
    complex(8) :: h
    complex(8), allocatable :: eps(:,:)
    complex(8), allocatable :: w2b(:), bw1(:)
    integer(4), allocatable :: ipiv(:)
    complex(8), allocatable :: work(:)
    
    character(len=10) :: sname="calcinveps"        

!!EXTERNAL ROUTINES: 
    external zgetrf, zgetri
    external zhetrf, zhetri
    complex(8), external :: zdotc, zdotu
!!REVISION HISTORY:

! Created 31.01.2007 by JH
! Revisited Dec 2011 by DIN
!
!EOP
!BOC
    call timesec(tstart)
    
    ! local arrays for body and wings
    allocate(eps(mbsiz,mbsiz))
    if (Gamma) allocate(bw1(mbsiz),w2b(mbsiz))

    ! LAPACK working arrays
    lwork = 64*mbsiz
    allocate(ipiv(mbsiz))
    allocate(work(lwork))

    ! lopp over frequencies
    do iom = iomstart, iomend
    
      ! array for body and its inverse
      eps(1:mbsiz,1:mbsiz) = epsilon(1:mbsiz,1:mbsiz,iom)

      select case (freq%fconv)
        case('refreq')
          call zgetrf(mbsiz,mbsiz,eps,mbsiz,ipiv,info)
          call errmsg0(info,sname,"calling zgetrf")
          call zgetri(mbsiz,eps,mbsiz,ipiv,work,lwork,info)
          call errmsg0(info,sname,"calling zgetri")
        case('imfreq')
          call zhetrf('u',mbsiz,eps,mbsiz,ipiv,work,lwork,info)
          call errmsg0(info,sname,"calling zhetrf")
          call zhetri('u',mbsiz,eps,mbsiz,ipiv,work,info)
          call errmsg0(info,sname,"calling zhetri")
      end select

!----------------------------------------------------------------------!
!                    Gamma point                                       !
!----------------------------------------------------------------------!

      if (Gamma) then
       
        ! Isotropic averaging along q0
        q0eps(:) = input%gw%scrcoul%q0eps(:)
        modq0    = q0eps(1)**2+q0eps(2)**2+q0eps(3)**2
        if (modq0 > 1.d-8) q0eps(:) = q0eps(:)/modq0

        h = epsh(iom,1,1)*q0eps(1) + &
        &   epsh(iom,2,2)*q0eps(2) + &
        &   epsh(iom,3,3)*q0eps(3)
        eps00(iom,1,1) = h

        epsw1(:,iom,1) = epsw1(:,iom,1)*q0eps(1) + &
        &                epsw1(:,iom,2)*q0eps(2) + &
        &                epsw1(:,iom,3)*q0eps(3)
        
        epsw2(:,iom,1) = epsw2(:,iom,1)*q0eps(1) + &
        &                epsw2(:,iom,2)*q0eps(2) + &
        &                epsw2(:,iom,3)*q0eps(3)

        select case (freq%fconv)
          case('refreq')
            call zgemv('n',mbsiz,mbsiz,zone,eps,mbsiz, &
            &           epsw1(:,iom,1),1,zzero,bw1,1)
            call zgemv('t',mbsiz,mbsiz,zone,eps,mbsiz, &
            &           epsw2(:,iom,1),1,zzero,w2b,1)
          case('imfreq')
            call zhemv('u',mbsiz,zone,eps,mbsiz, &
            &           epsw1(:,iom,1),1,zzero,bw1,1)
            call zhemv('u',mbsiz,zone,eps,mbsiz, &
            &           epsw2(:,iom,1),1,zzero,w2b,1)
            w2b = conjg(w2b)
        end select
 
        ! head^{-1}
        h = 1.d0 / (h - zdotu(mbsiz,epsw2(:,iom,1),1,bw1,1))
        epsh(iom,1,1)  = h

        ! wing^{-1}
        epsw1(:,iom,1) = -h*bw1(:)
        epsw2(:,iom,1) = -h*w2b(:)

        ! body^{-1}
        do jm = 1, mbsiz
        do im = 1, mbsiz
          eps(im,jm) = eps(im,jm) + epsw1(im,iom,1)*epsw2(jm,iom,1) / h
        enddo
        enddo

      endif ! q = 0

      !--------------------------------------------
      ! store into the same array
      !--------------------------------------------
      epsilon(1:mbsiz,1:mbsiz,iom) = eps(1:mbsiz,1:mbsiz)

    enddo ! iom
    
    deallocate(eps)
    if (Gamma) deallocate(bw1,w2b)
    deallocate(ipiv,work)
    
    !===========================================
    ! \epsilon^{-1}_{ij}-\delta_{ij}  
    !===========================================
    do iom = iomstart, iomend
      if (Gamma) epsh(iom,1,1) = epsh(iom,1,1)-zone
      do im = 1, mbsiz
        epsilon(im,im,iom) = epsilon(im,im,iom)-zone
      end do  
    end do ! iom
    
    call timesec(tend)
    time_dfinv = time_dfinv+tend-tstart

end subroutine
!EOC
