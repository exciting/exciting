!BOP
!
! !ROUTINE: calcvxcnn
!
! !INTERFACE:
      subroutine calcvxcnn

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $v^{xc}_{nn}(\vec{k})$, 
!only for valence states.
!
!
! !USES:
      use modinput
      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: ik
      integer(4) :: ist
      integer(4) :: ia, is
      integer(4) :: ngp
      
      complex(8) :: zt1
      complex(8), allocatable :: apwalm(:,:,:,:)
      complex(8), allocatable :: evecfv(:,:)
      complex(8), allocatable :: evec(:)
      complex(8), allocatable :: h(:)

      real(8)    :: tstart, tend
!
! !EXTERNAL ROUTINES: 
!
      complex(8), external :: zdotc
      
      external getevecfv
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
      allocate(vxcraa(apwordmax,0:input%groundstate%lmaxmat,apwordmax,0:input%groundstate%lmaxapw,0:lmmaxvr,natmtot))
      if (allocated(vxcrloa)) deallocate(vxcrloa)
      allocate(vxcrloa(nlomax,apwordmax,0:input%groundstate%lmaxmat,0:lmmaxvr,natmtot))
      if (allocated(vxcrlolo)) deallocate(vxcrlolo)
      allocate(vxcrlolo(nlomax,nlomax,0:lmmaxvr,natmtot))
!      
      call vxcrad
!      
      allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
      allocate(evecfv(nmatmax,nstfv))
      allocate(evec(nmatmax))
      if (allocated(vxcnn)) deallocate(vxcnn)
      allocate(vxcnn(nstfv,nkpt))

      do ik = 1, nkpt
        
        ngp=ngk(1,ik)

!       get the eigenvectors from file
        call getevecfvgw(idikp(ik),evecfv)

!       find the matching coefficients
        call match(ngp,gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)

        do ist = 1, nstfv

          allocate(h(nmat(1,ik)))
          h(:)=zzero

!         muffin-tin contributions
          do is = 1, nspecies
            do ia = 1, natoms(is)

              call vxcaa (.True., is, ia, ngp, apwalm, evecfv(:,ist), h)
              call vxcalo (.True., is, ia, ngp, apwalm, evecfv(:,ist), h)
              call vxclolo (.True., is, ia, ngp, evecfv(:,ist), h)

            end do
          end do

!         interstitial contributions
          call vxcistl(.True., ngp, igkig(:,1,ik), evecfv(:,ist), h)
        
          vxcnn(ist,ik)=zdotc(nmat(1,ik),evecfv(:,ist),1,h,1)
        
          deallocate(h)
        
        end do ! ist
        
      enddo ! ik  

!     print results into file VXCNN.OUT
      call writevxcnn

      deallocate(apwalm)
      deallocate(evecfv)
      deallocate(evec)
       
      deallocate(vxcraa)
      deallocate(vxcrloa)
      deallocate(vxcrlolo)
      
      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'CALCVXCNN')
          
      end subroutine calcvxcnn
!EOC      


