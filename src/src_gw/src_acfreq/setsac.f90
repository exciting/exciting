!BOP
!
! !ROUTINE: setsac
!
! !INTERFACE: 
      subroutine setsac(iopac,nomeg,npar,en,sc,omega,apar,poles)

! !DESCRIPTION:
!
! This subroutine set up the parameters for the selfenery analytic continuation (SAC) 
!
!
! !USES:

      implicit none
      integer(4),intent(in):: iopac         ! Option for which method to use 
                                            ! 0  -- Thiele's reciprocal difference method as described in
                                            !       H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977)
                                            ! 1  -- Rojas, Godby and Needs (PRL 74, 1827 (1996)

      integer(4),intent(in):: nomeg         ! Number of frequency points along the imaginary axis
      integer(4),intent(in):: npar          ! Number of parameters (== 2*npol (npol= the number of poles)  
      real(8),intent(in):: en               ! LDA eigen-energy around which an analyic continuation is performed
      real(8),intent(in):: omega(nomeg)     ! frequency points along the imaginary frequency 
      complex(8),intent(in) :: sc(nomeg)    ! correlation selfenergy along the imaginary frequency 
      complex(8),intent(out):: apar(npar)   ! fitted paramters to calculate selfenergy 
      complex(8),intent(out):: poles(npar)  ! the positions of poles of fitted selfenergy 
      
     
! !LOCAL VARIABLES:   
      integer(4) :: iw,step,ipar,ierr
      integer(4) :: iwpa(npar)
      real(8) :: varsq    ! Square root of the mean square error
      real(8)    :: xin(nomeg),anl(2*npar)      ! Values of the function
      complex(8) :: yin(nomeg)             ! Values of the function
      complex(8) :: cx(npar),cy(npar)        ! input for Pade's approximation 
      complex(8) :: coefs(0:npar/2)
      logical::lpolish


! !INTRINSIC ROUTINES: 
      
      intrinsic dble

! !EXTERNAL ROUTINES: 
      external setrgn
      external setpatrd

! !REVISION HISTORY:
! Created: Nov. 19, 2007 by  JH
!EOP
!BOC
       poles=0.d0

       if(en.gt.0.0d0)then
         xin(1:nomeg)=omega(1:nomeg)
         yin(1:nomeg)=sc(1:nomeg)
       else  
         xin(1:nomeg)=-omega(1:nomeg)
         yin(1:nomeg)=conjg(sc(1:nomeg))
       endif
!
! Choose the frequency points for the Pade's approximation 
!
       call setwpa(nomeg,npar,omega,iwpa)       

       do ipar=1,npar
         iw=iwpa(ipar)
         cx(ipar)=cmplx(0.d0,xin(iw))
         cy(ipar)=yin(iw)
       enddo

       if(iopac.eq.0 .or. npar .eq. nomeg)  then 
         call setpatrd(npar,cx,cy,apar)

       elseif(iopac.eq.2) then 
         do ipar=1,npar
           iw=iwpa(ipar)
           cx(ipar)=cmplx(-en,xin(iw))
           cy(ipar)=yin(iw)
         enddo
         call setpatrd(npar,cx,cy,apar)

       else 
         call init_c(cx,cy,apar,npar)

!         coefs(0)=cmplx(1.d0,0.d0)
!         coefs(1:npar/2)=apar(npar/2+1:npar)
!         call zroots(coefs,npar/2,poles,lpolish,ierr)
!         if(ierr.ne.0) then
!           write(6,*) "setsac: WARNING --- ZROOTS for poles are not &
!     & converged (0)"
!         endif

         do ipar=1,npar
           anl(ipar)=real(apar(ipar))
           anl(ipar+npar)=aimag(apar(ipar))
         enddo
         call nllsq(xin,yin,nomeg,anl,2*npar,varsq)
         do ipar=1,npar
           apar(ipar)=cmplx(anl(ipar),anl(ipar+npar),8)
         enddo

!         coefs(0)=cmplx(1.d0,0.d0)
!         coefs(1:npar/2)=apar(npar/2+1:npar)
!         call zroots(coefs,npar/2,poles(npar/2+1:npar),lpolish,ierr)
!         if(ierr.ne.0) then 
!           write(6,*) "setsac: WARNING --- ZROOTS for poles are not &
!     & converged (1) "
!         endif 
       endif 
       return
          
      endsubroutine 
!EOC        
