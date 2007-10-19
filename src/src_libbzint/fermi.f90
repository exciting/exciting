!BOP
!
! !ROUTINE: fermi
!
! !INTERFACE:
      subroutine fermi(nik,nbd,eband,ntet,tetc,wtet,vt,nel,sp,efer)
!
! !DESCRIPTION:
!  This subroutine calculated the Fermi energy with the tetrahedron method
!
! !USES:
      
      implicit none     
      
! !INPUT PARAMETERS:

      integer(4), intent(in) :: nik ! Number of irreducible k-points
      
      integer(4), intent(in) :: nbd ! Maximum number of bands
      
      real(8), intent(in) :: eband(nik,*) ! Band energies
      
      integer(4), intent(in) :: ntet        ! Number of tetrahedra
      
      integer(4), intent(in) :: tetc(4,*)! id. numbers of the corners
!                                             of the tetrahedra
  
      integer(4), intent(in) :: wtet(*)  ! weight of each tetrahedron
      
      real(8), intent(in)    :: vt         ! the volume of the tetrahedra

      real(8), intent(in)    :: nel        ! number of electrons
      
      logical, intent(in)    :: sp         ! .true. for spin polarized
!                                             case
      
! !OUTPUT PARAMETERS:      
      
      real(8), intent(out)   :: efer       ! the fermi energy

! !REVISION HISTORY:
!
! Created 10th. March 2004 by RGA
!
! !LOCAL VARIABLES:

      integer(4) :: ik, ib, it

      real(8) :: fact,emin,emax,eint,ocmin,ocmax,ocint,df
      
      real(8), external :: idos
      real(8), external :: dostet
      

!EOP
!BOC

      fact=1.0d0
      if(.not.sp)fact = fact + 1.0d0
      
      
      emin=100.00
      emax=-100.00
      do ib=1,nbd
        do ik=1,nik
!          write(25,*)eband(ik,ib)
          if(eband(ik,ib).lt.emin)emin=eband(ik,ib)
          if((eband(ik,ib).gt.emax).and.(eband(ik,ib).lt.5.0d+2))       &
   &                emax=eband(ik,ib)
        enddo
      enddo
      
      ocmin=fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,emin)
      ocmax=fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,emax)
      
      write(25,*)'emin =', emin,'ocmin =', ocmin
      write(25,*)'emax =', emax,'ocmax =', ocmax
           
      if(ocmax.le.nel) stop 'not enough bands'
      
      ocint=0.0d0
      it=0      
      do while((dabs(ocint-nel).ge.1.0d-4).and.((emax-emin).gt.1.0d-4)  &
     &    .and.(it.lt.100))
        it=it+1
        eint=emin+(emax-emin)/2
        ocint = fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,eint)
!        write(25,*)'emin emax eint'
        
        write(25,1)it,emin, emax, eint
!        write(25,*)'ocmin ocmax ocint'
        write(25,2)it,ocmin, ocmax, ocint
        if(ocint.gt.nel)then
          emax=eint
          ocmax=ocint
        else
          emin=eint
          ocmin=ocint
        endif
      enddo
      it=0      
      write(25,*)
      do while((dabs(ocint-nel).ge.1.0d-8).and.((emax-emin).gt.1.0d-8)  &
     &    .and.(it.lt.200))
        it=it+1
        eint=emin+(emax-emin)/(ocmax-ocmin)*(nel-ocmin)
        ocint = fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,eint)
!        write(25,*)'emin emax eint'
        write(25,1)it,emin, emax, eint
!        write(25,*)'ocmin ocmax ocint'
        write(25,2)it,ocmin, ocmax, ocint
        if(ocint.gt.nel)then
          emax=eint
          ocmax=ocint
        else
          emin=eint
          ocmin=ocint
        endif
      enddo
      df=dostet(nik,nbd,eband,ntet,tetc,wtet,vt,eint)
!      write(25,*)'df = ', df
      if(df.lt.1.0d-9)then
        emin=maxval(eband(1:nik,1:nbd),mask=eband(1:nik,1:nbd).lt. &
     &              eint)
        emax=minval(eband(1:nik,1:nbd),mask=eband(1:nik,1:nbd).gt. &
     &              eint)
        eint=0.5d0*(emin+emax)
      endif        
      efer=eint
      write(25,*)'Fermi energy = ',efer, df
    1 format('it =',i4,' emin =',f18.10,' emax =',f18.10,' eint =',f18.10) 
    2 format('it =',i4,' ocmin =',f18.10,' ocmax =',f18.10,' ocint =',f18.10) 
      end subroutine fermi
!EOC        
          

      
      
            



      
