
subroutine calcepsilon_ppm(iq,iomstart,iomend)

    use modinput
    use modmain
    use modgw
    implicit none
    ! input parameters
    integer, intent(in) :: iq
    integer(4), intent(in) :: iomstart, iomend
    ! local variables
    integer :: iom, im, jm
    real(8) :: freqs_
    complex(8) :: om, omp, tom, R, A
    complex(8) :: eh(2)
    complex(8), allocatable :: ew1(:,:), ew2(:,:)
    complex(8), allocatable :: eb(:,:,:)
    
    complex(8), parameter :: ieta=cmplx(0.d0,1.d-6,8)
    
    freqs_ = freq%freqs(2)
    ! first frequency is always \omega->0
    ! set the input plasmon frequency the second in the grid
    freq%freqs(2) = input%gw%scrcoul%omegap
    
    ! calculate the dielectric tensor at \omega=0 and \omega=\omega_p
    call init_dielectric_function(mbsiz,1,2,Gamma)
    call calcepsilon(iq,1,2)
    call calcinveps(1,2)
    
    ! save \epsilon in local arrays for PPM fitting
    allocate(eb(mbsiz,mbsiz,1:2))
    eb(:,:,1:2) = epsilon(:,:,1:2)
    if (Gamma) then
      eh(1:2) = epsh(1:2,1,1)
      allocate(ew1(mbsiz,1:2),ew2(mbsiz,1:2))
      ew1(:,1:2) = epsw1(:,1:2,1)
      ew2(:,1:2) = epsw2(:,1:2,1)
    end if
    
    ! restore back the original frequency grid
    freq%freqs(2) = freqs_
    call delete_dielectric_function(Gamma)
    
    ! calculate the dielectric tensor in PPM approximation for the rest
    call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
    omp = cmplx(input%gw%scrcoul%omegap**2,0.d0,8)
    do iom = iomstart, iomend
      om = cmplx(freq%freqs(iom)**2,0.d0,8)
      !-------
      ! Body
      !-------
      do jm = 1, mbsiz
        do im = 1, mbsiz
            ! \tilde{\omega}^2
            tom = eb(im,jm,2)/(eb(im,jm,1)-eb(im,jm,2))*omp
            ! \Omega^2
            A = -eb(im,jm,1)*tom
            ! PPM: \epsilon^{-1}-1
            epsilon(im,jm,iom) = -A/(om+tom) 
        end do
      end do
      if (Gamma) then
        !-------
        ! Head
        !-------
        ! \tilde{\omega}^2
        tom = eh(2)/(eh(1)-eh(2))*omp
        ! \Omega^2
        A = -eh(1)*tom
        ! PPM: head^{-1}-1
        epsh(iom,1,1) = -A/(om+tom)
        !-------
        ! Wings
        !-------
        do im = 1, mbsiz
          ! \tilde{\omega}^2
          tom = ew1(im,2)/(ew1(im,1)-ew1(im,2))*omp
          ! \Omega^2
          A = -ew1(im,1)*tom
          ! PPM: wing1^{-1}
          epsw1(im,iom,1) = -A/(om+tom)
          !-----------------
          ! \tilde{\omega}^2
          tom = ew2(im,2)/(ew2(im,1)-ew2(im,2))*omp
          ! \Omega^2
          A = -ew2(im,1)*tom
          ! PPM: wing2^{-1}
          epsw2(im,iom,1) = -A/(om+tom)
        end do
      end if
    end do ! iom
        
    deallocate(eb)
    if (Gamma) deallocate(ew1,ew2)
    
    return
end subroutine
