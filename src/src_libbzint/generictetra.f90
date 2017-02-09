!BOP
!
! !ROUTINE: generictetra
!
! !INTERFACE:

      subroutine generictetra(corners,wt_out,ical,info)
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

      use tetra_internal   , only : sgnfrq, omgga,fout,ztol_vol,&
     &               vol_small_tetra,iop_integ,ldbg_bzint,tol_taylor

      use polyhedron
      use order
    
! !INPUT PARAMETERS:

      implicit none
      
      real(8), intent(in) :: corners(3,4) ! Coordinates of the four nodes
      integer(4), intent(in) :: ical

! !OUTPUT PARAMETERS:            

      integer(4), intent(out) :: info
      real(8), intent(out) :: wt_out(4)      ! The four weights corresponding
!                                         to the original coordinates

! !LOCAL VARIABLES:

      integer(4) :: i,j,k

      integer(4), dimension(4) :: ind

      integer(4) :: sigeq

      real(8)    :: vol, det, max_de_small
      real(8)    :: pi

      real(8), dimension(3,3) :: vec 
      
      real(8), dimension(4) :: delta_e_big_tet, delta_e_small_tet,      &
     &                         wt_small_tet, wt_tmp
!
! !EXTERNAL ROUTINES: 
!
      external ksurf
      external sorteq
      external stweight_imag
      external stweight_itaylor
      external stweight_real
      external stweight_rtaylor
      
! !REVISION HISTORY:
!
! Created 6. april 2004 by RGA, 
! last revised by RGA on May.2008.
!
! Intrinsic functions

      det(i,j)=vec(2,i)*vec(3,j)-vec(2,j)*vec(3,i)
!EOP
!BOC
      pi=4.0d0*atan(1.0)
      info=0
      wt_out(1:4)=0.0d0
      wt_small_tet(1:4)=0.0d0
      wt_tmp(1:4)=0.0d0
!
! Calculate the energy differences
!
      do i=1,4
        delta_e_big_tet(i)=f(i)-e(i)
      enddo
!
! Calculate the volume of the small tetrahedron
!
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
      vol_small_tetra = vol 
      if(vol.lt.ztol_vol) goto 999
!
! For frequency dependent weights, calculate the energy diferences at the
! corners of the small tetrahedron and store the maximum absolute value
!
      if (sgnfrq.ne.1)then
        max_de_small=0.0d0
        do i=1,4
          delta_e_small_tet(i) = (delta_e_big_tet(2)-delta_e_big_tet(1))&
     &      * corners(1,i) + (delta_e_big_tet(3) - delta_e_big_tet(1))  &
     &      * corners(2,i) + (delta_e_big_tet(4) - delta_e_big_tet(1))  &
     &      * corners(3,i) + delta_e_big_tet(1)
          if(abs(delta_e_small_tet(i)).gt.max_de_small)max_de_small=abs(delta_e_small_tet(i))
        enddo
      endif  

      select case (sgnfrq)

      case(1)
        wt_small_tet(1:4)=1.0d0/2.40d+1
    
      case(2)
        if(iop_integ.eq.1) then 
          call stweight_numeric(0,delta_e_small_tet,omgga,wt_small_tet)
        else 
          if (omgga.gt.tol_taylor*max_de_small)then
            call stweight_rtaylor(delta_e_small_tet,omgga,wt_small_tet)
          else
            call sorteq(delta_e_small_tet,ind,sigeq)
            call stweight_real(delta_e_small_tet,omgga,sigeq,wt_tmp)
            do i=1,4
              wt_small_tet(ind(i))=wt_tmp(i)
            enddo
          endif
        endif 
      case(3)
        if(iop_integ.eq.1) then 
          call stweight_numeric(1,delta_e_small_tet,omgga,wt_small_tet)
          if(ldbg_bzint) write(fout,100) wt_small_tet 
        else 
          if (omgga.gt.tol_taylor*max_de_small)then
            call stweight_itaylor(delta_e_small_tet,omgga,wt_small_tet)
          else    
            call sorteq(delta_e_small_tet,ind,sigeq)
            call stweight_imag(delta_e_small_tet,omgga,sigeq,wt_tmp)
            do i=1,4
              wt_small_tet(ind(i))=wt_tmp(i)
            enddo
          endif 
        endif  
  100 format('#wt=',4g16.6)
 
      case(4)
        call sort(4,delta_e_small_tet,ind)
        call ksurf(delta_e_small_tet,omgga,wt_tmp)
        do i=1,4
          wt_small_tet(ind(i))=-1.0d0*pi*wt_tmp(i)
        enddo
      end select

      wt_out(1)=sum(wt_small_tet(1:4))
      do i=1,3
        do j=1,4
          wt_out(i+1)=wt_out(i+1)+wt_small_tet(j)*corners(i,j)
        enddo
        wt_out(1)=wt_out(1)-wt_out(i+1)
      enddo
      do i=1,4
        wt_out(i)=wt_out(i)*vol
      enddo  
      
      return
      
999   if(ical.ne.1)then ! If ical != 1 there is an error, write a warning
        info = 1
        write(fout,*)'WARNING:four vertices nearly in the same plane'
        write(fout,*)'ical = ',ical
        write(fout,*) 'vol = ',vol
        write(fout,*)'nodes'
        do i=1,4
          write(fout,*)corners(1:3,i)
        enddo  
        write(fout,*)'vecs'
        do i=1,3
          write(fout,*)vec(i,1:3)
        enddo 
      endif ! ical != 1
      wt_out(1:4)=0.0d0
!
!    If the frequency is zero, the contribution from the Fermi surface has
!    to be calculated         
!
      if((sgnfrq.eq.4).and.(omgga.lt.1.0d-12)) then
        call ksurf(e,ef,wt_small_tet)
        do i=1,4
          wt_out(i)=-1.0d0*wt_small_tet(i)*pi
        enddo
      endif
      
      end subroutine generictetra
!EOC      
