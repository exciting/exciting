!
subroutine calcinveps(iomstart,iomend)
!
! This subroutine calculates the inverse of the dielectric function
!
    use modmain
    use modgw
    use mod_mpi_gw, only : myrank

    implicit none
    integer(4), intent(in) :: iomstart, iomend
    integer(4) :: i, j, iom
    integer(4) :: im, jm
    integer(4) :: info, lwork
    real(8)    :: tstart, tend
    real(8)    :: wto, wlo
    complex(8) :: h, w, f
    complex(8), allocatable :: eps(:,:)
    complex(8), allocatable :: w2b(:), bw1(:)
    integer(4), allocatable :: ipiv(:)
    complex(8), allocatable :: work(:)

    character(len=10) :: sname="calcinveps"

    external zgetrf, zgetri
    external zhetrf, zhetri
    complex(8), external :: zdotc, zdotu

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
          call zpotrf('u', mbsiz, eps, mbsiz, info )
          call errmsg0(info,sname,"calling zpotrf")
          call zpotri('u', mbsiz, eps, mbsiz, info )
          call errmsg0(info,sname,"calling zpotri")
      end select

!----------------------------------------------------------------------!
!                    Gamma point                                       !
!----------------------------------------------------------------------!

      if (Gamma) then

        ! Isotropic averaging (alternative - spherical averaging, see XS/angavsc0)
        eps00(1,1,iom) = (epsh(1,1,iom)+epsh(2,2,iom)+epsh(3,3,iom)) / 3.0
        epsw1(:,1,iom) = (epsw1(:,1,iom)+epsw1(:,2,iom)+epsw1(:,3,iom)) / sqrt(3.0)
        epsw2(:,1,iom) = (epsw2(:,1,iom)+epsw2(:,2,iom)+epsw2(:,3,iom)) / sqrt(3.0)

        ! Block-matrix inversion
        call zgemv('n', mbsiz, mbsiz, zone, eps, mbsiz, &
        &           epsw1(:,1,iom), 1, zzero, bw1, 1)
        call zgemv('t', mbsiz, mbsiz, zone, eps, mbsiz, &
        &           epsw2(:,1,iom), 1, zzero, w2b, 1)

        ! head^{-1}
        h = 1.d0 / (eps00(1,1,iom) - zdotu(mbsiz, epsw2(:,1,iom), 1, bw1, 1))
        epsh(1,1,iom)  = h

        ! wing^{-1}
        epsw1(:,1,iom) = -h*bw1(:)
        epsw2(:,1,iom) = -h*w2b(:)

        ! body^{-1}
        do jm = 1, mbsiz
          do im = 1, mbsiz
            eps(im,jm) = eps(im,jm) + h*bw1(im)*w2b(jm)
          enddo
        enddo

        ! if (input%gw%eph == 'polar') then
        !   wlo = 2.99d-3
        !   wto = 1.39d-3
        !   w = zi*freq%freqs(iom)
        !   ! f = (wlo**2-wto**2)/(wlo**2-w**2)
        !   f = (wlo**2-wto**2) / wlo**2
        !   do jm = 1, mbsiz
        !     do im = 1, mbsiz
        !       eps(im,jm) = eps(im,jm) - f*epsw1(im,1,iom)*epsw2(jm,1,iom)/epsh(1,1,iom)
        !     end do
        !   end do
        !   epsw1(:,1,iom) = epsw1(:,1,iom) - epsw1(:,1,iom)*f
        !   epsw2(:,1,iom) = epsw2(:,1,iom) - epsw2(:,1,iom)*f
        ! end if

      endif ! q = 0

      !--------------------------------------------
      ! store into the same array
      !--------------------------------------------
      epsilon(1:mbsiz,1:mbsiz,iom) = eps(1:mbsiz,1:mbsiz)

    enddo ! iom

    deallocate(ipiv, work)
    deallocate(eps)
    if (Gamma) deallocate(bw1, w2b)

    !===========================================
    ! \epsilon^{-1}_{ij}-\delta_{ij}
    !===========================================
    do iom = iomstart, iomend
      if (Gamma) epsh(1,1,iom) = epsh(1,1,iom)-zone
      do im = 1, mbsiz
        epsilon(im,im,iom) = epsilon(im,im,iom)-zone
      end do
    end do ! iom

    call timesec(tend)
    time_dfinv = time_dfinv+tend-tstart

end subroutine
