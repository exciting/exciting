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

    complex(8), allocatable :: b(:,:)
    integer(4), allocatable :: ipiv(:)
    complex(8), allocatable :: work(:)
    
!!EXTERNAL ROUTINES: 
    external zgetrf, zgetri
    external zhetrf, zhetri
       
!!REVISION HISTORY:

! Created 31.01.2007 by JH
! Revisited Dec 2011 by DIN
!
!EOP
!BOC
    call timesec(tstart)
    
    ! local arrays for body and wings
    allocate(b(mbsiz,mbsiz))

    ! LAPACK working arrays
    lwork = 64*mbsiz
    allocate(ipiv(mbsiz))
    allocate(work(lwork))

    ! lopp over frequencies
    do iom = iomstart, iomend
    
      ! array for body and its inverse
      b(1:mbsiz,1:mbsiz) = epsilon(1:mbsiz,1:mbsiz,iom)

      !------------------------
      ! invert the body
      !------------------------
      ! ATTENTION: zhetrf works (input/output) only with UPPER matrix triangle 
      !call zhetrf('u',mbsiz,b,mbsiz,ipiv,work,lwork,info)
      !call errmsg0(info,'calcinveps',"calling zhetrf")
      !call zhetri('u',mbsiz,b,mbsiz,ipiv,work,info)
      !call errmsg0(info,'calcinveps',"calling zhetri")
      call zgetrf(mbsiz,mbsiz,b,mbsiz,ipiv,info)
      call errmsg0(info,'calcinveps','calling zgetrf')
      call zgetri(mbsiz,b,mbsiz,ipiv,work,lwork,info)
      call errmsg0(info,'calcinveps','calling zgetri')
      
      ! \Gamma point treatment
      if (Gamma) call angular_averaging(iom,mbsiz,b)
      
      !--------------------------------------------
      ! store into the same array
      !--------------------------------------------
      epsilon(1:mbsiz,1:mbsiz,iom) = b(1:mbsiz,1:mbsiz)

    enddo ! iom
    
    deallocate(b)
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
