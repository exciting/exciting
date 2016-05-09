!BOP
!
!!ROUTINE: setbarcev
!
!!INTERFACE:
!
subroutine setbarcev(evtol)
!
! !DESCRIPTION:
!
! This subroutine reset the eigenvectors and eigenvalues of 
! bare coulomb matrix in terms of evtol
!
!!USES:
    use modmain,               only : zone, zzero
    use modgw
    use mod_mpi_gw,            only : myrank

!!INPUT PARAMETERS: 
    implicit none
    real(8), intent(in) :: evtol 

!!LOCAL VARIABLES:
    integer(4) :: im, jm
    integer(4) :: immax
    real(8) :: test1, test2
    complex(8) :: vc
    integer(4), allocatable :: im_kept(:) ! indicate which barc eigenvectors
                                          ! are kept as basis functions
    complex(8), allocatable :: wi0new(:)

!!REVISION HISTORY:
!
! Created July 31,2009 by Hong Jiang
! Readjusted Jan, 2012 by DIN
!
!EOP
!BOC

    ! Reduce the basis size by choosing eigenvectors of barc with 
    ! eigenvalues larger than evtol
    allocate(im_kept(matsiz))
    im_kept(:) = 1
    mbsiz = matsiz
    do im = 1, matsiz
      if (barcev(im) < evtol) then 
        im_kept(im) = 0
        mbsiz = mbsiz-1
      end if 
    end do 

    if (Gamma) then
      call calcwmix0
      allocate(wi0new(matsiz))
      call zgemv('c',matsiz,matsiz,zone,vmat,matsiz,wi0,1,zzero,wi0new,1)
      ! find the index of the diagonalized barc eigenvector that has maximal 
      ! overlap with G=0 (constant) plane wave
      test2 = 0.0d0
      do im = 1, matsiz
        test1 = dble(wi0new(im)*conjg(wi0new(im)))
        if (test1 > test2) then
          immax = im
          test2 = test1
        end if
      end do
      if (myrank==0) then
        write(fgw,*)'- Maximum singular eigenvector ###'
        write(fgw,10) immax, test2, barcev(immax)
      end if    
      ! exclude this matrix element      
      if (im_kept(immax) == 1) then 
        im_kept(immax) = 0
        mbsiz = mbsiz-1
      end if
      deallocate(wi0,wi0new)
    end if
    10 format("immax, max(wi0new), barcev(immax): ",i4,4x,f8.3,4x,f9.6)
    
    if (mbsiz < matsiz) then
      if (myrank==0) then
        write(fgw,*) "Info(setbarcev): Product basis size has been changed"
        write(fgw,*) " - Old basis set size =", matsiz
        write(fgw,*) " - New basis set size =", mbsiz
      end if
    end if
    
    ! Build the trasformation matrix
    if (allocated(barc)) deallocate(barc)
    allocate(barc(matsiz,mbsiz))
    barc(:,:) = zzero
    
    im = 0
    do jm = 1, matsiz 
      if (im_kept(jm) == 1) then 
        im = im+1
        vc = cmplx(barcev(jm),0.d0,8)
        barc(:,im) = vmat(:,jm)*sqrt(vc)
      end if 
    end do   
    deallocate(im_kept)

    return
end subroutine
!EOC
