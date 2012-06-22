
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: generictetra
!
! !INTERFACE:

      subroutine generictetra(corners,ww2,ical,info)
!
! !DESCRIPTION:
!
! This subroutine calculates the integrals:
! 
! For the case of normal weight $sigfreq=1$:
! \begin{equation}
! \begin{align}
! w(1)=&\iiint\limits_T (1-x-y-z) dx dy dz \\
! w(2)=&\iiint\limits_T x dx dy dz \\
! w(3)=&\iiint\limits_T y dx dy dz \\ 
! w(4)=&\iiint\limits_T z dx dy dz 
! \end{align}
! \end{equation}
! 
! For the case of the weight including real frequency contribution, 
! $sigfreq=2$:
!
! \begin{eqnarray}
! w(1)=\iiint\limits_T \frac{1-x-y-z}{\omega-\epsilon_{n'k-q}+\epsilon_{nk}} 
! dx dy dz \nonumber \\
! w(2)=\iiint\limits_T \frac{x}{\omega-\epsilon_{n'k-q}+\epsilon_{nk}} 
! dx dy dz \nonumber \\
! w(3)=\iiint\limits_T \frac{y}{\omega-\epsilon_{n'k-q}+\epsilon_{nk}} 
! dx dy dz \nonumber \\
! w(4)=\iiint\limits_T \frac{z}{\omega-\epsilon_{n'k-q}+\epsilon_{nk}} 
! dx dy dz
! \end{eqnarray}
!
! The whole contribution is got by setting $\omega$ to $-\omega$ and call this
! program again.
!
! For the case of weight including imaginary frequency, $sigfreq=3$:
!
! \begin{eqnarray}
! w(1)=\iiint\limits_T \frac{-2(1-x-y-z)\omega}{\omega^2+(\epsilon_{n'k-q}+\epsilon_{nk})^2}
! dx dy dz \nonumber \\
! w(2)=\iiint\limits_T \frac{-2x\omega}{\omega^2+(\epsilon_{n'k-q}+\epsilon_{nk})^2} 
! dx dy dz \nonumber \\
! w(3)=\iiint\limits_T \frac{-2y\omega}{\omega^2+(\epsilon_{n'k-q}+\epsilon_{nk})^2} 
! dx dy dz \nonumber \\
! w(4)=\iiint\limits_T \frac{-2z\omega}{\omega^2+(\epsilon_{n'k-q}+\epsilon_{nk})^2} 
! dx dy dz
! \end{eqnarray}
!
! where $T$ is a tetrahedron of corners \verb"nodes(i,j)".

! For the case of surface integration weight, $sigfreq=4$:
!
! \begin{eqnarray}
! w(1)=\iiint\limits_T (1-x-y-z)\Delta(\epsilon_{nk}-\epsilon_{n'k-q}+\omega)
! dx dy dz \nonumber \\
! w(2)=\iiint\limits_T x\Delta(\epsilon_{nk}-\epsilon_{n'k-q}+\omega)
! dx dy dz \nonumber \\
! w(3)=\iiint\limits_T y\Delta(\epsilon_{nk}-\epsilon_{n'k-q}+\omega)
! dx dy dz \nonumber \\
! w(4)=\iiint\limits_T z\Delta(\epsilon_{nk}-\epsilon_{n'k-q}+\omega)
! dx dy dz
! \end{eqnarray}
!
! where $T$ is a tetrahedron of corners \verb"nodes(i,j)".



! !USES:

      use tetra_internal   , only : sgnfrq, omgga
!<sag>
      use control, only : restype, tetradbglv
!</sag>

      use polyhedron

      use order
    
! !INPUT PARAMETERS:

      implicit none
      
      real(8), intent(in) :: corners(3,4) ! Coordinates of the four nodes
      integer(4), intent(in) :: ical

! !OUTPUT PARAMETERS:            

      integer(4), intent(out) :: info
      real(8), intent(out) :: ww2(4)      ! The four weights corresponding
!                                         to the original coordinates

! !LOCAL VARIABLES:

      integer(4) :: i,j,k,sigeq

      integer(4), dimension(4) :: ind

      integer(4), dimension(5) :: sig
     
      real(8)    :: vol, det,  weit1, weit2, maxv

      real(8), dimension(3,3) :: vec 
      
      real(8), dimension(4) :: energydif, v, ww1, ww0, w1

      ! <sag>
      ! resonant/antiresonant parts separation, 2007/05/01 by S. Sagmeister
      logical :: tres,tares
      ! </sag>
!
! !EXTERNAL ROUTINES: 
!
      external ksurf
 
! !REVISION HISTORY:
!
! Created 6. april 2004 by RGA, last revised by XZL on Dec.7th.
!
! Intrinsic functions

      !<sag>
      ! the line below for "det" is an intrinsic function definition (?)
      !</sag>
      det(i,j)=vec(2,i)*vec(3,j)-vec(2,j)*vec(3,i)
!EOP
!BOC
      !<sag>
      ! zero the variables "weit1" and "weit2" because of separate treatment
      ! of resonant and antiresonant weights
      weit1=0.d0
      weit2=0.d0
      ! assign resonance type
      tres=(restype.eq.0).or.(restype.eq.1)
      tares=(restype.eq.0).or.(restype.eq.2)
      !</sag>

      info=0
      ww2(1:4)=0.0d0

      do i=1,4
        energydif(i)=f(i)-e(i)
      enddo
      
      v(1:4)=0.0d0

      select case (sgnfrq)

      case(1)
        do i=1,3
          do j=1,3
            vec(i,j)=corners(j,i+1)-corners(j,1)
          enddo
        enddo
        vol=0.0d0
        do i=1,3
          j=mod(i,3)+1
          k=mod(j,3)+1
          vol=vol+vec(1,i)*det(j,k)
        enddo
        vol=abs(vol)/6.0d0
        if(vol.eq.0.0)goto 999
        ww2(1)=vol
        do i=1,3
          do j=1,4
            ww2(i+1)=ww2(i+1)+2.50d-1*vol*corners(i,j)
          enddo
          ww2(1)=ww2(1)-ww2(i+1)
        enddo
    
      case(2)

        do i=1,4
          v(i)=(energydif(2)-energydif(1))*corners(1,i)
          v(i)=v(i)+(energydif(3)-energydif(1))*corners(2,i)
          v(i)=v(i)+(energydif(4)-energydif(1))*corners(3,i)
          v(i)=v(i)+energydif(1)
        enddo


!------------------------------------------------------------
!       Added by RGA 2.12.04 Begin
!------------------------------------------------------------
        do i=1,3
          do j=1,3
            vec(i,j)=corners(j,i+1)-corners(j,1)
          enddo
        enddo
        vol=0.0d0
        do i=1,3
          j=mod(i,3)+1
          k=mod(j,3)+1
          vol=vol+vec(1,i)*det(j,k)
        enddo
        vol=abs(vol)
        if(vol.eq.0.0)goto 999
!------------------------------------------------------------
!       Added by RGA 2.12.04 End
!------------------------------------------------------------
        maxv=0.0d0
        do i=1,4
          if(abs(v(i)).gt.maxv)maxv=abs(v(i))
        enddo

        if (omgga.ge.1.0d+1*maxv)then
          !<sag>
          !<commented>
!!$          call redifwtaylor(v,omgga,ww1)
!!$          call redifwtaylor(v,-omgga,ww2)
          !</commented>
          if (tres) call redifwtaylor(v,omgga,ww1)
          if (tares) call redifwtaylor(v,-omgga,ww2)
          !</sag>
          do i=1,4
           ww1(i)=ww1(i)+ww2(i)
          enddo
        else
          call sorteq(v,ind,sig)

          sigeq=sig(5)

          ww2(1:4)=0.0d0

          ww1(1:4)=0.0d0

          w1(1:4)=0.0d0

          !<sag>
          !<commented>
!!$          call redifwt(v,omgga,sigeq,weit1)
!!$
!!$          call redifwt(v,-omgga,sigeq,weit2)
          !</commented>
          if (tres) call redifwt(v,omgga,sigeq,weit1)
          if (tares) call redifwt(v,-omgga,sigeq,weit2)
          !</sag>

          ww1(ind(1))=weit1+weit2
 
          !<sag>
          !<commented>
!!$          call redifwtz(v,omgga,sigeq,weit1)
!!$
!!$          call redifwtz(v,-omgga,sigeq,weit2) 
          !</commented>
          if (tres) call redifwtz(v,omgga,sigeq,weit1)
          if (tares) call redifwtz(v,-omgga,sigeq,weit2) 
          !</sag>

          ww1(ind(4))=weit1+weit2

          !<sag>
          !<commented>
!!$          call redifwtx(v,omgga,sigeq,weit1)
!!$
!!$          call redifwtx(v,-omgga,sigeq,weit2)
          !</commented>
          if (tres) call redifwtx(v,omgga,sigeq,weit1)
          if (tares) call redifwtx(v,-omgga,sigeq,weit2)
          !</sag>

          ww1(ind(2))=weit1+weit2

          !<sag>
          !<commented>
!!$          call redifwty(v,omgga,sigeq,weit1)
!!$ 
!!$          call redifwty(v,-omgga,sigeq,weit2)
          !</commented>
          if (tres) call redifwty(v,omgga,sigeq,weit1)
          if (tares) call redifwty(v,-omgga,sigeq,weit2)
          !</sag>

          ww1(ind(3))=weit1+weit2 
        endif
        ww2(:)=0.0d0
        ww2(1)=ww1(1)+ww1(2)+ww1(3)+ww1(4)
        do i=1,3
          do j=1,4
            ww2(i+1)=ww2(i+1)+ww1(j)*corners(i,j)
          enddo
          ww2(1)=ww2(1)-ww2(i+1)
        enddo
!RGA ADED 2.12.04 begin
        do i=1,4
          ww2(i)=ww2(i)*vol
        enddo  
!RGA ADED 2.12.04 end

      case(3)
        
        do i=1,4
          v(i)=(energydif(2)-energydif(1))*corners(1,i)
          v(i)=v(i)+(energydif(3)-energydif(1))*corners(2,i)
          v(i)=v(i)+(energydif(4)-energydif(1))*corners(3,i)
          v(i)=v(i)+energydif(1)
        enddo


!------------------------------------------------------------
!       Added by RGA 2.12.04 Begin
!------------------------------------------------------------
        do i=1,3
          do j=1,3
            vec(i,j)=corners(j,i+1)-corners(j,1)
          enddo
        enddo
        vol=0.0d0
        do i=1,3
          j=mod(i,3)+1
          k=mod(j,3)+1
          vol=vol+vec(1,i)*det(j,k)
        enddo
        vol=abs(vol)
        if(vol.eq.0.0)goto 999
!------------------------------------------------------------
!       Added by RGA 2.12.04 End
!------------------------------------------------------------
        maxv=0.0d0
        do i=1,4
          if(abs(v(i)).lt.1.0d-10)v(i)=0.0d0
          if(abs(v(i)).gt.maxv)maxv=abs(v(i))
        enddo
        
        if (omgga.ge.1.0d+1*maxv)then
          call edifwtaylor(v,omgga,ww1)
        else    
          call sorteq(v,ind,sig)

          sigeq=sig(5)

          ww2(1:4)=0.0d0
 
          ww1(1:4)=0.0d0

          w1(1:4)=0.0d0

          call edifwt(v,omgga,sigeq,weit1)

          ww1(ind(1))=weit1
 
          call edifwtz(v,omgga,sigeq,weit1)
 
          ww1(ind(4))=weit1

          call edifwtx(v,omgga,sigeq,weit1)

          ww1(ind(2))=weit1

          call edifwty(v,omgga,sigeq,weit1)
 
          ww1(ind(3))=weit1
        
        endif  
 
        ww2(1)=ww1(1)+ww1(2)+ww1(3)+ww1(4)


        do i=1,3
          do j=1,4
            ww2(i+1)=ww2(i+1)+ww1(j)*corners(i,j)
          enddo
          ww2(1)=ww2(1)-ww2(i+1)
        enddo
!RGA ADED 2.12.04 begin
        do i=1,4
          ww2(i)=ww2(i)*vol
        enddo  
!RGA ADED 2.12.04 end




!--------------------------------------------------------------------
!     case 4:    Added by XZL for surface integration 12.05.05 Begin
!--------------------------------------------------------------------        

      case(4)

        do i=1,4
          v(i)=(energydif(2)-energydif(1))*corners(1,i)
          v(i)=v(i)+(energydif(3)-energydif(1))*corners(2,i)
          v(i)=v(i)+(energydif(4)-energydif(1))*corners(3,i)
          v(i)=v(i)+energydif(1)
        enddo
    
        call sort(4,v,ind)
        do i=1,3
          do j=1,3
            vec(i,j)=corners(j,i+1)-corners(j,1)
          enddo
        enddo
        vol=0.0d0
        do i=1,3
          j=mod(i,3)+1
          k=mod(j,3)+1
          vol=vol+vec(1,i)*det(j,k)
        enddo
        vol=abs(vol)
        if(vol.lt.1.0d-15)goto 999

        call ksurf(v,omgga,ww0)

        do i=1,4
          ww1(ind(i))=ww0(i)
        enddo
        ww2(1)=ww0(1)+ww0(2)+ww0(3)+ww0(4)
        ww2(2:4)=0.0d0
        do i=1,3
          do j=1,4
            ww2(i+1)=ww2(i+1)+ww1(j)*corners(i,j)
          enddo
          ww2(1)=ww2(1)-ww2(i+1)
        enddo

        do i=1,4
          ww2(i)=ww2(i)*vol*(0.0d0-3.1415926d0)
        enddo
!--------------------------------------------------------------------
!     case 4:    Added by XZL for surface integration 12.05.05 End
!--------------------------------------------------------------------



      end select


      
      return
      
999   info = 1
!<sag>
      if (tetradbglv.gt.10) then
!</sag>
      write(*,*)'the four points are in the same plane'
      write(*,*)'ical = ',ical
      write(*,*)'nodes'
      do i=1,4
        write(*,*)corners(1:3,i)
      enddo  
      write(*,*)'vecs'
      do i=1,3
        write(*,*)vec(i,1:3)
      enddo  
!<sag>
      end if
!</sag>
      do i=1,4
       ww2(i)=0.0d0
      enddo
      if((sgnfrq.eq.4).and.(omgga.lt.1.0d-12)) then
        call ksurf(e,ef,ww2)
        do i=1,4
          ww2(i)=ww2(i)*(0.0d0-3.1415926d0)
        enddo
      endif
      
      end subroutine generictetra
!EOC      
