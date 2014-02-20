!BOP
!
! !MODULE: q0barc
      module q0barc

! !DESCRIPTION:
!
! This module declares the variables used in barcq0 and calcsing
!
! !PUBLIC VARIABLES:

        integer(4) :: nglen
        real(8),    allocatable :: glen0(:)
        real(8),    allocatable :: sing(:,:,:)
        complex(8), allocatable :: phase(:,:,:)
        
      end module q0barc
!EOP               
!---------------------------------------------------------------
!BOP
!
! !ROUTINE: barcq0
!
! !INTERFACE:
      subroutine barcq0

! !DESCRIPTION:
!
! This subroutine calculates the matrix of the bare coulomb potential for
! q=0 and atomic functions with L=0
!
! !USES:

      use modmain
      use modgw
      use modmpi
      use q0barc
!
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: is        ! Indexes the atoms (inequivalent) 
!                               for the columns of barc
      integer(4) :: ias       ! Indexes all the atoms 
!                               for the columns of barc
      integer(4) :: ia        ! Indexes the equivalent atoms 
!                               for the columns of barc
      integer(4) :: igl
                              ! imix(ipw) = locmatsiz + iipw)
      integer(4) :: imix      ! Indexes the mixed wave function 
!                               for the columns of barc (barc(imix,jmix))
      integer(4) :: irm       ! Indexes the mixed basis functions 
!                               of atom iat
      integer(4) :: js        ! Indexes the atoms (inequivalent) 
!                               for the files of barc
      integer(4) :: jas       ! Indexes all the atoms
!                               for the files of barc
      integer(4) :: ja        ! Indexes the equivalent atoms 
!                               for the files of barc
      integer(4) :: jmix      ! Indexes the mixed wave function 
!                               for the files of barc (barc(imix,jmix)).
      integer(4) :: jrm       ! Indexes the mixed basis functions 
!                               of atom js
      integer(4) :: l1        ! Angular momentum quantum number of the 
!                               atomic function irm
      integer(4) :: l2        ! Angular momentum quantum number of the 
!                               atomic function jrm
      integer(4) :: m1        ! z-component angular momentum quantum 
!                               number of the atomic function irm
      integer(4) :: m2        ! z-component angular momentum quantum 
!                               number of the atomic function jrm
      
      real(8) :: qg1len       ! |qg1|
      real(8) :: invglen2
      real(8) :: len4
      real(8) :: t1, t2, t3
 
      complex(8) ::  sumipw
 
!
! !REVISION HISTORY:
! 
! Created 16th. March 2004 by RGA
! Last modified 31. March 2005 by RGA
! Revisited: June 2011 by DIN
!
!EOP
!BOC
      call cpu_time(t1)
      call calcsing
      call cpu_time(t2)

      imix=0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          do irm=1,nmix(ias)
            l1=bigl(ias,irm)
            if(l1.eq.0)then
              imix=imix+1
              jmix=0
              do js=1,nspecies
                do ja=1,natoms(js)
                  jas=idxas(ja,js)
                  do jrm = 1, nmix(jas)
                    l2=bigl(jas,jrm)
                    if(l2.eq.0)then
                      jmix=jmix+1
                      sumipw=zzero
                      do igl = 2, nglen
                        qg1len = glen0(igl)
                        len4=qg1len*qg1len*qg1len*qg1len
                        invglen2=1.0d0/len4
                        sumipw=sumipw+phase(ias,jas,igl)*                &
     &                         cmplx(sing(ias,irm,igl)*sing(jas,jrm,igl)*&
     &                         invglen2,0.0d0,8)
                      enddo ! igl  
                      barc(imix,jmix)=16.0d0*pi*pi*vi*sumipw
                    else  
                      do m2=-l2,l2
                        jmix=jmix+1
                      enddo ! m2
                    endif  
                  enddo ! jrm
                enddo ! ja
              enddo ! js
            else  
              do m1=-l1,l1
                imix = imix + 1
              enddo ! m1
            endif    
          enddo ! irm
        enddo ! ia
      enddo ! is                
      deallocate(sing)
      deallocate(phase)
      deallocate(glen0)

      call cpu_time(t3)
      if(t1.lt.0.0d0)write(fgw,*)'warning, t1 < 0'
      if(t2.lt.0.0d0)write(fgw,*)'warning, t2 < 0'
      if(t3.lt.0.0d0)write(fgw,*)'warning, t3 < 0'
      if (rank == 0) then
        call write_cputime(fgw,t3-t1,'    BARCQ0')
        call write_cputime(fgw,t2-t1,'    CALCSING')
      end if

      return
      end subroutine barcq0
!EOC      
!BOP
!
! !ROUTINE: calcsing
!
! !INTERFACE:
      subroutine calcsing
!
! !DESCRIPTION: 
!
!
!
! !USES:
!
      use modinput
      use modmain
      use modgw
      use q0barc
	use modmpi
! !LOCAL VARIABLES:
!
      implicit none

      integer(4) :: is
      integer(4) :: ias
      integer(4) :: ia    
      integer(4) :: ng1,ng2,ng3
      integer(4) :: igl
      integer(4) :: ig1,ig2,ig3
      integer(4) :: ipw
      integer(4) :: irm ! (Counter) runs over radial mixed basis
!                         functions
      integer(4) :: irp ! (Counter), runs over the radial mesh points.
      integer(4) :: js
      integer(4) :: jas
      integer(4) :: ja
      integer(4) :: jgl
      integer(4) :: nrp    ! Number of radial (logarithmic) 
!                                        mesh points for atom is
      integer(4) :: l   ! Angular momentum quantum number of the 
!                         mixed atomic function irm.
      integer(4) :: ng

      integer(4), dimension(3) :: igvec ! Integer coordinates of the G-vector
      integer(4), allocatable :: gind(:,:),gind4(:)
      
      real(8) :: qglen ! Length of q+G
      real(8) :: gmax
      real(8) :: kk
      real(8) :: glprev
      real(8) :: gpr
      real(8) :: x

      real(8), dimension(3) :: gvec ! Coordinates of the vector G
      real(8), dimension(3) :: dpos
      real(8), allocatable :: fr(:), gr(:) ! Temporary arrays for the radial
!                                         functions
      real(8), allocatable :: cf(:,:)
      real(8), allocatable :: sinf(:)    
      real(8), allocatable :: rp(:)     ! Radial mesh points
      real(8), allocatable :: glen(:)     ! Radial mesh points
      
      complex(8) :: expg

! !INTRINSIC ROUTINES: 

      intrinsic dsqrt

!
! !REVISION HISTORY:
!
! Created: 17th. March 2004 by MF
! Last Modified: 30th. March 2004 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC
!
      gmax=10*input%gw%MixBasis%gmb*gkmax
      ng1=idint(gmax*pia(1))+1
      ng2=idint(gmax*pia(2))+1
      ng3=idint(gmax*pia(3))+1
      ng = (2*ng1+1)*(2*ng2+1)*(2*ng3+1)
      
      allocate(glen(ng))
      allocate(gind(3,ng))
      allocate(gind4(ng))

      ipw=0
      do ig1=-ng1,ng1
        igvec(1)=ig1
        do ig2=-ng2,ng2
          igvec(2)=ig2
          do ig3=-ng3,ng3
            igvec(3)=ig3
!           Transform igvec to cartesian coordinates
            gvec=matmul(bvec,dble(igvec))
!           Calculate the length of the qpg:
            kk=sqrt(gvec(1)*gvec(1)+gvec(2)*gvec(2)+gvec(3)*gvec(3))
            if(kk.lt.gmax)then
              ipw=ipw+1
              glen(ipw)=kk
              gind(1:3,ipw)=igvec(1:3)
            endif
          enddo
        enddo
      enddo

      ng = ipw
!     sort by increasing length using shell algorithm
      call shelsort(ng,gind(:,1:ng),glen(1:ng))
      
      glprev=-1.0d0
      igl=0
      do ipw=1,ng
        if(abs(glen(ipw)-glprev).gt.1.0d-10)then
          igl=igl+1
          glprev=glen(ipw)
        endif
        gind4(ipw)=igl
      enddo
      nglen=igl
      
      allocate(glen0(nglen))
      glen0=0.d0
      jgl=0
      do ipw=1,ng
        if(gind4(ipw).ne.jgl)then
          jgl=gind4(ipw)
          glen0(jgl)=glen(ipw)
        endif
      enddo    
      deallocate(glen)
!     
      allocate(phase(natmtot,natmtot,nglen))
      phase=zzero
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias=idxas(ia,is)
          do js=1,nspecies
            do ja=1,natoms(is)
              jas=idxas(ja,js)
              dpos(1:3)= atposl(1:3,ia,is)-atposl(1:3,ja,js)
              do ipw = 2, ng
                igl=gind4(ipw)
                igvec(1:3)=gind(1:3,ipw)
                gpr=dble(igvec(1))*dpos(1)+dble(igvec(2))*dpos(2)+ &
               &    dble(igvec(3))*dpos(3)
                expg=cmplx(cos(2.0d0*pi*gpr),-sin(2.0d0*pi*gpr),8)
                phase(ias,jas,igl)=phase(ias,jas,igl)+expg
              enddo ! ipw
            enddo ! ja
          enddo ! js
        enddo ! ia
      enddo ! is
      deallocate(gind,gind4)

      allocate(sing(natmtot,maxnmix,nglen)) 
      sing=0.0d0
      
      do is=1,nspecies
        nrp=nrmt(is)
        allocate(rp(nrp))          
        allocate(sinf(nrp))          
        allocate(fr(nrp))
        allocate(gr(nrp))
        allocate(cf(3,nrp))
        rp(1:nrp)=spr(1:nrp,is)
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          do igl=2,nglen
            qglen=glen0(igl)
!
!           Calculate the spherical bessel function at each mesh point
!
            do irp = 1, nrp
              x=rp(irp)*qglen
              sinf(irp)=sin(x)
            enddo   ! irp
            do irm=1,nmix(ias)
              l=bigl(ias,irm)
              if(l.eq.0)then
                do irp=1,nrp
                  fr(irp)=umix(ias,irm,irp)*sinf(irp)
                enddo  
!               Integrate the wavefunctions:
                call fderiv(-1,nrp,rp,fr,gr,cf)
                sing(ias,irm,igl) = gr(nrp)
              endif  
            enddo ! irm
          enddo ! igl
        enddo ! ia  

        deallocate(rp)
        deallocate(fr)
        deallocate(gr)
        deallocate(cf)
        deallocate(sinf)
      enddo !is
      
      end subroutine calcsing
!EOC
