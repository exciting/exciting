!BOP
! 
! !ROUTINE: calcrlamint
!
! !INTERFACE:
      subroutine calcrlamint()

! !DESCRIPTION:
!
! This subroutine calculates the set of integrals:
!
!\begin{equation}
!\verb"rtl(ias,i)"\equiv \left\langle r^bl\right\rangle_{aNL} = %
!\int\limits_0^{R^a_{MT}}(r^a)^{bl+2}\upsilon_{aNL}(r^a)dr^a
!\end{equation}
!
!and 
!
!\begin{equation}
!\verb"rint(ias,N,N',bl)"\equiv \left\langle%
!\begin{array}{c}
!  r_<^{bl} \\
!  r_>^{bl +1} \\
!\end{array}%
!\right\rangle_{aNL,N'bl}
!=\iint\limits_0^{R^a_{MT}}%
!\upsilon_{aNL}(r^a_1)\frac{r_<^bl}{r_>^{bl+1}}\upsilon_{aNL}(r_2^a)%
!(r_1^a)^2dr_1^a(r_2^a)^2dr_2^a
!\end{equation}
!
!that will be needed for the calculation of the bare coulomb potential
!
!
! !USES:

      use modmain
      use modgw


! !LOCAL VARIABLES:

      implicit none

      integer(4) :: ia
      integer(4) :: ias  ! number of the corresponding 
!                          inequivalent atom.
      integer(4) :: irm  ! Counter: Runs over radial mixed basis functions.
      integer(4) :: jrm  ! Counter: Runs over radial mixed basis functions.
      integer(4) :: ijrm ! Joint index for (irm,jrm) pairs for compressed
!                          storage
      integer(4) :: is
      integer(4) :: bl   ! Angular momentum quantum number of the mixed
!                          basis functions
      integer(4) :: l1   ! just a counters
      integer(4) :: ir   ! Counter: Runs over radial mesh point.
!                          jrj(ias)

      real(8) :: rxov    ! value of the integral

      real(8), allocatable :: rrtol(:)
      real(8), allocatable :: ui(:)
      real(8), allocatable :: uj(:)
      real(8), allocatable :: fr(:)
      real(8), allocatable :: gr(:)
      real(8), allocatable :: cf(:,:)
!
! !REVISION HISTORY:
! 
! Created Dic. 2003 by RGA
! Last Modification: July. 20th. 2004 by RGA
! Revisited 5.05.2011 by DIN 
!
!EOP
!BOC
!     Initializations
      if (allocated(rtl)) deallocate(rtl)
      allocate(rtl(natmtot,maxnmix))
      if (allocated(rrint)) deallocate(rrint)
      allocate(rrint(natmtot,maxnmix*(maxnmix+1)/2))

      rtl=0.0d0
      rrint=0.0d0

!     Local arrays
      allocate(rrtol(nrmtmax))
      allocate(ui(nrmtmax),uj(nrmtmax))
      allocate(fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax))      

      do is=1,nspecies
        if(debug)write(701,100)is,trim(spname(is))
        do ia=1,natoms(is)
          ias=idxas(ia,is)

!-----------------------------------------------------------------
!                      Calculate rtl
!----------------------------------------------------------------- 
          if(debug)write(701,*)
          if(debug)write(701,101)

          rrtol(1:nrmt(is))=spr(1:nrmt(is),is)
          bl=0
          do irm=1,nmix(ias)
!
!           if bl has changed, calculate r^(bl+1) for each grid point
!
            if(bigl(ias,irm).gt.bl)then

              do l1=bl+1,bigl(ias,irm)
                do ir=1,nrmt(is)
                  rrtol(ir)=rrtol(ir)*spr(ir,is)
                enddo ! ir
              enddo ! l1

              bl=bigl(ias,irm)

            else if (bigl(ias,irm).lt.bl) then

              write(fgw,*)'WARNING!!!!'
              write(fgw,*)'radial mixed functions not ordered by ', &
     &                  'increasing bl'
              do l1=bigl(ias,irm),bl
                do ir=1,nrmt(is)
                  rrtol(ir)=rrtol(ir)/spr(ir,is)
                enddo
              enddo
              bl=bigl(ias,irm)  

            endif

            do ir=1,nrmt(is)
              fr(ir)=umix(ias,irm,ir)*rrtol(ir)
            enddo
            call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
            
            if (gr(nrmt(is)).lt.0.0d0) then
              rtl(ias,irm)=-1.0d0*gr(nrmt(is))
              umix(ias,irm,1:nrmt(is))=-1.0d0*umix(ias,irm,1:nrmt(is))
            else
              rtl(ias,irm)=gr(nrmt(is))
            endif
            
            if(debug)write(701,102)irm,rtl(ias,irm)

          enddo ! irm

!-----------------------------------------------------------------
!                      Calculate rrint 
!----------------------------------------------------------------- 
          rrtol(1:nrmt(is))=spr(1:nrmt(is),is)
          bl=0
          ijrm=0
          if(debug)write(701,20)

          do irm=1,nmix(ias)
!
!           if bl has changed, calculate r^(bl+1) for each grid point 
!
            if (bigl(ias,irm).gt.bl) then
              
              do l1=bl+1,bigl(ias,irm)
                do ir=1,nrmt(is)
                  rrtol(ir)=rrtol(ir)*spr(ir,is)
                enddo ! ir
              enddo ! l1
              bl=bigl(ias,irm)
            
            else if (bigl(ias,irm).lt.bl) then
              write(fgw,*)'WARNING!!!!'
              write(fgw,*)'radial mixed functions not ordered by ', &
             &          'increasind bl'
              do l1=bigl(ias,irm),bl
                do ir=1,nrmt(is)
                  rrtol(ir)=rrtol(ir)/spr(ir,is)
                enddo ! ir
              enddo ! l1
              bl=bigl(ias,irm)  
            
            endif
!
!           store the corresponding mixed function in a local array
!
            ui(1:nrmt(is))=umix(ias,irm,1:nrmt(is))
            uj(1:nrmt(is))=ui(1:nrmt(is))
            ijrm=irm+(irm*(irm-1))/2
            
            call drinteg(is,rxov)
            rrint(ias,ijrm)=rxov
            
            if(debug)write(701,21)irm,irm,bigl(ias,irm), &
           &                      bigl(ias,irm),rrint(ias,ijrm)
!          
!           now loop over jrm for the double integrals
!
            do jrm=irm+1,nmix(ias)
              ijrm=irm+(jrm*(jrm-1))/2
              if(bigl(ias,irm).eq.bigl(ias,jrm))then
                uj(1:nrmt(is))=umix(ias,jrm,1:nrmt(is))
                call drinteg(is,rxov)
              else
                rxov=0.0d0
              endif
              rrint(ias,ijrm)=rxov
              
              if(debug)write(701,21)irm,jrm,bigl(ias,irm),bigl(ias,jrm),      &
             &                      rrint(ias,ijrm)
            enddo ! jrm

          enddo ! irm

        enddo ! ia
      enddo ! is
      deallocate(ui,uj,rrtol,fr,gr,cf)   
      
   20 format(/,5x,'Double radial integrals',/,12x,'N1',2x,'N2',2x,'L1',&
     &       2x,'L2',5x,'rrint')
   21 format(10x,4i4,g18.10)
  100 format("Species : ",I4,", ",A)
  101 format(5x,'Radial integrals',/,13x,'N',3x,'<r^L|v_L>',9x,        &
     &      '<r^(L+2)|v_L> ')
  102 format(10x,i4,(1pg18.10,1x))
      return

      CONTAINS

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      subroutine drinteg(is,xint)
!
! The procedure used consists in writing it as:
!
!\begin{equation}
!\begin{aligned}
!&\int\limits_{0}^{R_{MT}}\upsilon_{aNL}(r_1)%
!\Biggl[r_1^{-(L+1)}\int\limits_{0}^{r_1}{r_2^{L}%
!\upsilon_{aN'L}(r_2)\left(r_2\right)^2dr_2}-\\
!&-r_1^{L}\int\limits_{R_{MT}}^{r_1}{{r_2}^{-L-1}%
!\upsilon_{aN'L}(r_2)\left(r_2\right)^2dr_2}\Biggr]
!\left(r_1\right)^2dr_1
!\end{aligned}
!\end{equation}
!
! Using an integration subroutine that return the intermediate values of
! the integral, reduces the calculation to just three integrations
!
! !SEE ALSO: fint.f
!
!

      implicit none
      integer(4), intent(in) :: is  ! Inequivalent atom for which the integral is calculated

! !OUTPUT PARAMETERS:
!
      real(8), intent(out):: xint   ! value of the integral
!
! !DEFINED PARAMETERS:

      integer(4), parameter :: nsym = 6 ! Order of the polynom used to fit
!                                         the function and make the
!                                         analytical integration

!
! !LOCAL VARIABLES:
!
      integer(4) :: irp ! runs over spherical grid points
      
      integer(4) :: jri                 ! grid size
      real(8), allocatable :: rr(:)     ! grid point array
      
      real(8), allocatable :: rtlint(:) ! 1-d integral <r^l>, the vector contains the intermediate values, that is rtl(irp)=int_0^{r(irp)}u1 r^l u2 dr
      real(8), allocatable :: irtl1(:)  ! the one dimensionan integral
!                                         <r^(-l-1)>, the vector contains the
!                                         intermediate values, that is
!                                         rtl(irp)=int_RMT^{r(irp)}u1 r^(-l-1)
!                                                  u2 dr
      real(8), allocatable :: rr2(:)    ! the grid points in inverse order
      real(8), allocatable :: uitrl(:)  ! u_NL x r^(L+1)
      real(8), allocatable :: ujtrl(:)  ! u_N'L x r^(L+1)
      real(8), allocatable :: uiorl(:)  ! u_NL / r^L
      real(8), allocatable :: ujorl(:)  ! u_N'L / r^L
      real(8), allocatable :: uij(:)


! !EXTERNAL ROUTINES:

      external fint
!
! !REVISION HISTORY:
!
! Created Jan. 2004
! Last Modified: Feb. 2nd. 2004
!
!EOP
!BOC
!
!     Allocate the needed arrays
!
      jri=nrmt(is)
      allocate(rtlint(jri))
      allocate(irtl1(jri))
      allocate(rr(jri))
      allocate(rr2(jri))
      allocate(uitrl(jri))
      allocate(ujtrl(jri))
      allocate(uiorl(jri))
      allocate(ujorl(jri))
      allocate(uij(jri))
     
      rr(:)=spr(1:jri,is)
      do irp=1,jri
        rr2(irp)=rr(jri-irp+1)
        uitrl(irp)=ui(irp)*rrtol(irp)
        ujtrl(irp)=uj(irp)*rrtol(irp)
        uiorl(jri-irp+1)=ui(irp)*rr(irp)/rrtol(irp)
        ujorl(irp)=uj(irp)*rr(irp)/rrtol(irp)
      enddo
!
!     Integrate in r_2:
!
      call fint(nsym,jri,rr,uitrl,rtlint)
      call fint(nsym,jri,rr2,uiorl,irtl1)
      do irp=1,jri
        uij(irp)= ujorl(irp)*rtlint(irp)-ujtrl(irp)*irtl1(jri-irp+1)
      enddo
!
!     Integrate in r_1:
!
      call fint(nsym,jri,rr,uij,rtlint)
      xint=rtlint(jri)

      deallocate(rtlint)
      deallocate(irtl1)
      deallocate(rr2)
      deallocate(uitrl)
      deallocate(ujtrl)
      deallocate(uiorl)
      deallocate(ujorl)
      deallocate(uij)

      end subroutine drinteg

      end subroutine calcrlamint
!EOC

