!BOP
!
! !ROUTINE: setbarcev
!
! !INTERFACE:
      subroutine setbarcev(evtol)

! !DESCRIPTION:
!
! This subroutine reset the eigenvectors and eigenvalues of 
! bare coulomb matrix in terms of evtol
!
! !USES:

      use modmain
      use modgw
!
! !INPUT PARAMETERS: 

      implicit none
      real(8),    intent(in) :: evtol 
!
! !LOCAL VARIABLES:

      integer(4) :: im,jm      ! Indexes the mixed wave function for the columns of barc (barc(im,jm))
      integer(4) :: immax
      real(8)    :: test1,test2

      integer(4), allocatable :: im_kept(:) ! indicate which barc eigenvectors are kept as basis functions 
      complex(8), allocatable :: wi0new(:)  ! eigenvalues of sqrt(barc)
!
! !REVISION HISTORY:
!
! Created July 31,2009 by Hong Jiang
! Adapted Jan, 2012 by DIN
!
!EOP
!BOC
!
!    Reduce the basis size by choosing eigenvectors of barc with 
!    eigenvalues larger than evtol
!
      allocate(im_kept(matsiz),wi0new(matsiz))
      im_kept = 1
      mbsiz=matsiz
      do im=1,matsiz
        if( barcev(im).lt.evtol ) then 
          im_kept(im)=0
          mbsiz=mbsiz-1
        endif 
      enddo 
!
!     Construct new basis set for q=0 by removing
!     all eigenvectors with eigenvalues smaller than evtol  
!
      if (Gamma) then
        call calcwmix0
        call zgemv('c',matsiz,matsiz,zone,vmat,matsiz,wi0,1,zzero,wi0new,1)
        !! find the index of the diagonalized barc eigenvector that has maximal overlapw with
        !! G=0 (constant) plane wave
        test2=0.0d0
        do im=1,matsiz
          test1= real(wi0new(im)*conjg(wi0new(im)))
          if(test1.gt.test2)then
            immax=im
            test2=test1
          endif
        enddo
        if(evtol.gt.0)then
          write(fgw,*)'- Maximum singular eigenvector ###'
          write(fgw,100) immax,test2,barcev(immax)
        end if
        if(im_kept(immax).eq.1) then 
          im_kept(immax)=0
          mbsiz=mbsiz-1
        endif 
      endif
  100 format("immax,max(wi0new),barcev(immax)=",i4,f8.3,e10.3)

      
      if(allocated(barcvm))deallocate(barcvm)
      allocate(barcvm(matsiz,mbsiz))
      if(allocated(vbas))deallocate(vbas)
      allocate(vbas(matsiz,mbsiz))

      im=0
      do jm=1,matsiz 
        if(im_kept(jm).eq.1) then 
          im=im+1
          vbas(:,im)=vmat(:,jm)
          barcvm(:,im)=vmat(:,jm)*sqrt(barcev(jm))
        endif 
      enddo   

      deallocate(im_kept)

      if(evtol.gt.0)then
        write(fgw,*) " Reduce the basis size by choosing eigenvectors of barc"
        write(fgw,*) " with eigenvalues larger than evtol"
        write(fgw,*) " evtol =", evtol
        write(fgw,*) "  - Old basis set size =", matsiz
        write(fgw,*) "  - New basis set size =", mbsiz
      end if

      return
      end subroutine 
