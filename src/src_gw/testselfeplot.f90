!BOP
!
! !ROUTINE: testselfeplot
!
! !INTERFACE: 
      subroutine testselfeplot

! !DESCRIPTION:
!
! This subroutine plots the selfenergy as a function of the frequency
!
!
! !USES:
     
      use modmain
      use modgw

! !LOCAL VARIABLES:

      integer(4) :: ie
      integer(4) :: i
      integer(4) :: iom
      integer(4) :: ipar
      integer(4) :: nom
      integer(4) :: npar
      integer(4) :: ikp
      integer(4) :: iemin, iemax
      integer(4) :: ir

      real(8), allocatable :: fr(:)
      real(8)    :: enk
      
      complex(8), allocatable :: comega(:)
      complex(8), allocatable :: cy1(:),cy2(:),cy3(:)
      complex(8), allocatable :: dy(:)
      complex(8), allocatable :: selfe1(:,:),selfe2(:,:),selfe3(:,:)
      complex(8), allocatable :: selfed(:,:)
      complex(8), allocatable :: trash(:)    ! temporary array (unused)
      complex(8), allocatable :: trash2(:,:) ! temporary array (unused)
      complex(8), allocatable :: acpar(:)
      complex(8), allocatable :: a(:),poles(:) ! Parameters for AC     
!EOP
!BOC

!
!     Plot parameters
!
      ikp=input%gw%iik
      iemin=input%gw%ibmin
      iemax=input%gw%ibmax
      ir=1

!     work only for AC with the multipole fitting (MPF) scheme
      iopac=1
!      
!     Read the selfenergy from SELFC.OUT
!
      allocate(selfec(nstfv,nkpt,nomeg))      
      call readselfc
!
!     Output files
!      
      open(unit=71,file='realselfc.out',form='formatted',status='unknown')
      write(71,*)'### Real part of the correlation selfenergy for different types of AC'
      write(71,*)'### omega[eV]    ratfun     RGN    Pade'
      
      open(unit=72,file='imagselfc.out',form='formatted',status='unknown')
      write(72,*)'### Imaginary part of the correlation selfenergy for different types of AC'
      write(72,*)'### omega[eV]    ratfun     RGN    Pade'

!
!     Calculate the number of parameters of the analitic function for the selfenergy
!      
      npar=2*npol
      nom=2*nomeg+1
      
      allocate(comega(1:nom))
      allocate(fr(1:nom))
      allocate(cy1(1:nom))
      allocate(cy2(1:nom))
      allocate(cy3(1:nom))
      allocate(dy(1:nom))
      allocate(selfe1(iemin:iemax,1:nom))
      allocate(selfe2(iemin:iemax,1:nom))
      allocate(selfe3(iemin:iemax,1:nom))
      allocate(selfed(iemin:iemax,1:nom))
      allocate(trash(1:2*npar))
      allocate(trash2(1:2*npar,1:2*npar))
      allocate(a(2*npar))
      allocate(acpar(npar),poles(npar))

      comega(nomeg+1)=zzero  
      fr(nomeg+1)=0.0d0

      do iom=1,nomeg
        fr(iom+nomeg+1)=freqs(iom)
        fr(nomeg+1-iom)=-1.0d0*freqs(iom)
      enddo
      
      if(ir.eq.0)then
        do iom=1,nom
          comega(iom)=cmplx(fr(iom),0.0d0,8)
        enddo
      else  
        do iom=1,nom
          comega(iom)=cmplx(0.0d0,fr(iom),8)
        enddo
      endif

      do ie=iemin,iemax
        
        enk=evaldft(ie,ikp)
        
        ! RGN ac scheme
        call setsac(1,nomeg,npar,enk,selfec(ie,ikp,1:nomeg),freqs,acpar,poles)
        
        do ipar=1,npar
          a(ipar)=real(acpar(ipar))
          a(ipar+npar)=aimag(acpar(ipar))
        enddo
        
        ! perform AC
        do iom=1,nom
          call ratfun(comega(iom),a,cy1(iom),trash,trash2,npol-1)
          call getsac(1,nomeg,npar,enk,comega(iom),freqs,acpar,cy2(iom),dy(iom))
        enddo

        ! PADE ac scheme
        call setsac(0,nomeg,npar,enk,selfec(ie,ikp,1:nomeg),freqs,acpar,poles)
        do iom=1,nom
          call getsac(0,nomeg,npar,enk,comega(iom),freqs,acpar,cy3(iom),dy(iom))
        enddo

        selfe1(ie,1:nom)=cy1(1:nom)  
        selfe2(ie,1:nom)=cy2(1:nom)  
        selfe3(ie,1:nom)=cy3(1:nom)
        selfed(ie,1:nom)=dy(1:nom)  
  
      enddo

      do i=iemin,iemax
        do iom=1,nom
          write(71,10)fr(iom)*hev,real(selfe1(i,iom)),                 &
     &                real(selfe2(i,iom)),real(selfe3(i,iom))
          write(72,10)fr(iom)*hev,aimag(selfe1(i,iom)),                 &
     &                aimag(selfe2(i,iom)),aimag(selfe3(i,iom))
        enddo
        write(71,*)
        write(72,*)
      enddo    

      deallocate(comega)
      deallocate(fr)
      deallocate(cy1)
      deallocate(cy2)
      deallocate(cy3)
      deallocate(dy)
      deallocate(a)
      deallocate(acpar)
      deallocate(poles)
      deallocate(selfe1)
      deallocate(selfe2)
      deallocate(selfe3)
      deallocate(selfed)
      deallocate(trash)
      deallocate(trash2)
       
      close(71)  
      close(72)              

   10 format(4(' ',g18.10))
      return
      
      return
      
      end subroutine testselfeplot
!EOC      
