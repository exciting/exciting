!BOP
!
! !ROUTINE: orthog_corewf
!
! !INTERFACE:
      subroutine orthog_corewf

! !USES:

      use modmain
      use modgw

!      use constants, only: cfein1, cfein2
!      use core,      only: ncore,ocmax,symbl,lcore
!      use lapwlo,    only: loor,lomax
!      use radwf,     only: ucore,uscore,u,us,udot,usdot,ulo,uslo
!      use struk,     only: nrpt

!
! !INPUT PARAMETERS:
      
      implicit none 

!
! !DESCRIPTION:
!
! Orthogonalizes the core wave-functions to the valence ones.
!
! !LOCAL VARIABLES:

      integer(4) :: is,ia,ias
      integer(4) :: icore ! Counter: run over core states 
      integer(4) :: lc    ! Angular momentum quantum number of the core wf
      integer(4) :: npt   ! number of points in the radial mesh
      integer(4) :: io,ir,ilo

      real(8) :: normcore ! norm of the core wave function
            
      real(8), allocatable :: uc(:)
      real(8), allocatable :: uval(:)

      real(8) :: fr(spnrmax)
      real(8) :: gr(spnrmax) 
      real(8) :: cf(3,spnrmax)

! !EXTERNAL ROUTINES: 

      external fderiv      
    
!
!EOP
!BOC      

      do is=1,nspecies

        npt=nrmt(is)
        allocate(uc(npt),uval(npt))
        
        do ia=1,natoms(is)
          ias=idxas(ia,is)
      
          do icore = 1, ncore(is) 
            lc=spl(icore,is)
!
!           The norm of the core wave function ucore 
!
            uc(1:npt)=ucore(1:npt,1,icore,ias)
            do ir = 1, npt
               fr(ir) = uc(ir)*uc(ir)*spr(ir,is)**2.d0
            end do !ir
            call fderiv(-1,npt,spr(1:npt,is),fr,gr,cf)
            normcore=gr(npt)
            write(6,*) 'ias, icore, lc, normcore: ', ias, icore, lc, normcore
!
!           Orthogonalize to apwfr
!
            do io = 1, apword(lc,is)
              if (apwdm(io,lc,is).eq.0) then
                uval(1:npt) = apwfr(1:npt,1,io,lc,ias)
                call orthogonalize
              end if
            end do !io
!
!           Orthogonalize to local orbital
!        
            do ilo = 1, nlorb(is)
              if (lc.eq.lorbl(ilo,is)) then
                uval(1:npt) = lofr(1:npt,1,ilo,ias)
                call orthogonalize
              endif
            end do ! ilo
!
!           Restore core wave function
!
            ucore(1:npt,1,icore,ias) = uc(1:npt)
           
          enddo !icore
        end do !ia
      
        deallocate(uc,uval)  
      
      end do !is

      CONTAINS
      
        subroutine orthogonalize
      
          implicit none

          integer(4) :: irp     ! index of the radial mesh

          real(8) :: overlap_factor 
          real(8) :: renorm_factor
          real(8) :: overlap
          real(8) :: normafter
          real(8) :: normval
!
!         The norm of the valence radial function 
!
          do irp = 1, npt
             fr(irp) = uval(irp)*uval(irp)*spr(irp,is)**2.d0
          end do !irp
          call fderiv(-1,npt,spr(1:npt,is),fr,gr,cf)
          normval=gr(npt)
!
!         Overlap between valence and core
!          
          do irp = 1, npt
             fr(irp) = uval(irp)*uc(irp)*spr(irp,is)**2.d0
          end do !irp
          call fderiv(-1,npt,spr(1:npt,is),fr,gr,cf)
          overlap=gr(npt)
!
!         Orthogonalization (Gramm-Schmith like)
!         
          overlap_factor = overlap/normval
          do irp = 1, npt
            uc(irp) = uc(irp) - overlap_factor * uval(irp)
          enddo ! irp
!
!         Norm of the orthogonalized core function
!
          do irp = 1, npt
             fr(irp) = uc(irp)*uc(irp)*spr(irp,is)**2.d0
          end do !irp
          call fderiv(-1,npt,spr(1:npt,is),fr,gr,cf)
          normafter=gr(npt)
!
!         Renormalization
!          
          renorm_factor = sqrt(normcore/normafter)
          do irp = 1, npt
            uc(irp) = renorm_factor * uc(irp)
          enddo  

        end subroutine orthogonalize  
      
      end subroutine orthog_corewf
!EOC  
