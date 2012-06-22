!BOP
!
! !ROUTINE: calcvxcnn
!
! !INTERFACE:
      subroutine calcvxcnn

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $v^{xc}_{nn}(\vec{k})$ 
!of equation \ref{vxc-10}, only for valence states.
!
!
! !USES:
      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: ik
      integer(4) :: ist
      integer(4) :: np
      
      complex(8) :: zt1
      complex(8), allocatable :: apwalm(:,:,:,:)
      complex(8), allocatable :: evecfv(:,:)
      complex(8), allocatable :: evec1(:)
      complex(8), allocatable :: evec2(:)
      complex(8), allocatable :: vxcapw(:)
      
      real(8)    :: tstart, tend
!
! !EXTERNAL ROUTINES: 
!
      complex(8), external :: zdotc
      
      external getevecfv
      external genvxcapw
      external match      

! !INTRINSIC ROUTINES: 

      intrinsic cpu_time
     
! !REVISION HISTORY:
! 
! Created  8th. Aug. 2006 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!
!BOC

      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'

!     allocate exchange-correlation integral arrays
      if (allocated(vxcraa)) deallocate(vxcraa)
      allocate(vxcraa(maxapword,0:input%groundstate%lmaxmat,maxapword,0:input%groundstate%lmaxapw,0:lmmaxvr,natmtot))
      if (allocated(vxcrloa)) deallocate(vxcrloa)
      allocate(vxcrloa(nlomax,maxapword,0:input%groundstate%lmaxmat,0:lmmaxvr,natmtot))
      if (allocated(vxcrlolo)) deallocate(vxcrlolo)
      allocate(vxcrlolo(nlomax,nlomax,0:lmmaxvr,natmtot))
      
      call vxcrad
      
      allocate(apwalm(ngkmax,maxapword,lmmaxapw,natmtot))
      allocate(evecfv(nmatmax,nstfv))
      allocate(evec1(nmatmax))
      allocate(evec2(nmatmax))

      allocate(vxcnn(nstfv,nkpt))

      do ik = 1, nkpt
        
        np=npmat(1,ik)
        allocate(vxcapw(np))

!       get the eigenvectors from file
        call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)

!       find the matching coefficients
        call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik), &
       &     sfacgk(:,:,1,ik),apwalm)

        call genvxcapw(ik,apwalm,vxcapw)
        
        do ist=1,nstfv
          evec1(:)=evecfv(:,ist)
          call zhpmv('U',nmat(1,ik),zone,vxcapw,evec1,1,zzero,evec2,1)
          zt1=zdotc(nmat(1,ik),evec1,1,evec2,1)
          vxcnn(ist,ik)=zt1
        enddo !ist  
        
        deallocate(vxcapw)
      
      enddo ! ik  

      call writevxcnn

      deallocate(apwalm)
      deallocate(evecfv)
      deallocate(evec1)
      deallocate(evec2)
      
      deallocate(vxcraa)
      deallocate(vxcrloa)
      deallocate(vxcrlolo)
      
      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'CALCVXCNN')
          
      end subroutine calcvxcnn
!EOC      


