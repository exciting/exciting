
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: surfnodes
!
! !INTERFACE:
      subroutine surfnodes
      
! !DESCRIPTION:
!
!This subroutine calculates the coordinates of the intersections of three
!planes between the four planes that form the surface of the tetrahedron
!and the two planes that approximate the fermi surface within it for the
!bands at $\vec{k}$ and $\vec{k}-\vec{q}$, in internal coordinates of the
!tetrahedron $\vec{x}=(\xi,\eta,\zeta)$. The equations defining the planes are:
!
!\begin{subequations}
!\begin{align}
!\xi=&0\\
!\eta=&0\\
!\zeta=&0\\
!\xi+\eta+\zeta=&1\\
!(\varepsilon_2-\varepsilon_1)\xi+(\varepsilon_3-\varepsilon_1)\eta+%
!(\varepsilon_4-\varepsilon_1)\zeta=&\varepsilon_F-\varepsilon_1\\
!(\varepsilon'_2-\varepsilon'_1)\xi+(\varepsilon'_3-\varepsilon'_1)\eta+%
!(\varepsilon'_4-\varepsilon'_1)\zeta=&\varepsilon_F-\varepsilon'_1
!\end{align}
!\end{subequations}
!
! where $\varepsilon_i=\varepsilon_n(\vec{k}_i)$ and 
! $\varepsilon'_i=\varepsilon_n(\vec{k}_i-\vec{q})$.
!
! The expressions for the nodes, in the order they are calculated, which
!corresponds to increasing ntype ordering, are:
!
!\begin{subequations}
!\begin{align}
! \vec{x}_1=&(0,0,0)\\
! \vec{x}_2=&(0,0,1)\\
! \vec{x}_3=&(0,1,0)\\
! \vec{x}_4=&(1,0,0)\\
! \vec{x}_5=&(0,0,\tfrac{\varepsilon_F-\varepsilon_1}%
! {\varepsilon_4-\varepsilon_1})\\
! \vec{x}_6=&(0,\tfrac{\varepsilon_F-\varepsilon_1}%
! {\varepsilon_3-\varepsilon_1},0)\\
! \vec{x}_7=&(\tfrac{\varepsilon_F-\varepsilon_1}%
! {\varepsilon_2-\varepsilon_1},0,0)\\
! \vec{x}_8=&(0,\tfrac{\varepsilon_4-\varepsilon_F}%
! {\varepsilon_4-\varepsilon_3},\tfrac{\varepsilon_3-\varepsilon_F}%
! {\varepsilon_3-\varepsilon_4})\\
! \vec{x}_9=&(\tfrac{\varepsilon_4-\varepsilon_F}%
! {\varepsilon_4-\varepsilon_2},0,\tfrac{\varepsilon_2-\varepsilon_F}%
! {\varepsilon_2-\varepsilon_4})\\
! \vec{x}_{10}=&(\tfrac{\varepsilon_3-\varepsilon_F}%
! {\varepsilon_3-\varepsilon_2},\tfrac{\varepsilon_2-\varepsilon_F}%
! {\varepsilon_2-\varepsilon_3},0)\\
! \vec{x}_{11}=&(0,0,\tfrac{\varepsilon_F-\varepsilon'_1}%
! {\varepsilon'_4-\varepsilon'_1})\\
! \vec{x}_{12}=&(0,\tfrac{\varepsilon_F-\varepsilon'_1}%
! {\varepsilon'_3-\varepsilon'_1},0)\\
! \vec{x}_{13}=&(\tfrac{\varepsilon_F-\varepsilon'_1}%
! {\varepsilon'_2-\varepsilon'_1},0,0)\\
! \vec{x}_{14}=&(0,\tfrac{\varepsilon'_4-\varepsilon_F}%
! {\varepsilon'_4-\varepsilon'_3},\tfrac{\varepsilon'_3-\varepsilon_F}%
! {\varepsilon'_3-\varepsilon'_4})\\
! \vec{x}_{15}=&(\tfrac{\varepsilon'_4-\varepsilon_F}%
! {\varepsilon'_4-\varepsilon'_2},0,\tfrac{\varepsilon'_2-\varepsilon_F}%
! {\varepsilon'_2-\varepsilon'_4})\\
! \vec{x}_{16}=&(\tfrac{\varepsilon'_3-\varepsilon_F}%
! {\varepsilon'_3-\varepsilon'_2},\tfrac{\varepsilon'_2-\varepsilon_F}%
! {\varepsilon'_2-\varepsilon'_3},0)\\
! \vec{x}_{17}=&(0,%
!\tfrac{(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_4)-%
!(\varepsilon_F-\varepsilon_4)(\varepsilon_F-\varepsilon'_1)}%
!{(\varepsilon_4-\varepsilon_3)(\varepsilon'_4-\varepsilon'_1)-%
!(\varepsilon_4-\varepsilon_1)(\varepsilon'_4-\varepsilon'_3)},
!\tfrac{(\varepsilon_F-\varepsilon_3)(\varepsilon_F-\varepsilon'_1)-%
!(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_3)}%
!{(\varepsilon_4-\varepsilon_3)(\varepsilon'_4-\varepsilon'_1)-%
!(\varepsilon_4-\varepsilon_1)(\varepsilon'_4-\varepsilon'_3)})\\
! \vec{x}_{18}=&(%
!\tfrac{(\varepsilon_F-\varepsilon_4)(\varepsilon_F-\varepsilon'_1)-%
!(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_4)}%
!{(\varepsilon_2-\varepsilon_4)(\varepsilon'_2-\varepsilon'_1)-%
!(\varepsilon_2-\varepsilon_1)(\varepsilon'_2-\varepsilon'_4)},0,%
!\tfrac{(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_2)-%
!(\varepsilon_F-\varepsilon_2)(\varepsilon_F-\varepsilon'_1)}%
!{(\varepsilon_2-\varepsilon_4)(\varepsilon'_2-\varepsilon'_1)-%
!(\varepsilon_2-\varepsilon_1)(\varepsilon'_2-\varepsilon'_4)})\\
! \vec{x}_{19}=&(%
!\tfrac{(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_3)-%
!(\varepsilon_F-\varepsilon_3)(\varepsilon_F-\varepsilon'_1)}%
!{(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_1)-%
!(\varepsilon_3-\varepsilon_1)(\varepsilon'_3-\varepsilon'_2)},%
!\tfrac{(\varepsilon_F-\varepsilon_2)(\varepsilon_F-\varepsilon'_1)-%
!(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_2)}%
!{(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_1)-%
!(\varepsilon_3-\varepsilon_1)(\varepsilon'_3-\varepsilon'_2)},0)\\
! \vec{x}_{20}=&(%
!\tfrac{(\varepsilon_F-\varepsilon_3)(\varepsilon_F-\varepsilon'_4)-%
!(\varepsilon_F-\varepsilon_4)(\varepsilon_F-\varepsilon'_3)}%
!{(\varepsilon_3-\varepsilon_4)(\varepsilon'_3-\varepsilon'_2)-%
!(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_4)},
!\tfrac{(\varepsilon_F-\varepsilon_4)(\varepsilon_F-\varepsilon'_2)-%
!(\varepsilon_F-\varepsilon_2)(\varepsilon_F-\varepsilon'_4)}%
!{(\varepsilon_3-\varepsilon_4)(\varepsilon'_3-\varepsilon'_2)-%
!(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_4)},\nonumber\\
!&\tfrac{(\varepsilon_F-\varepsilon_2)(\varepsilon_F-\varepsilon'_3)-%
!(\varepsilon_F-\varepsilon_3)(\varepsilon_F-\varepsilon'_2)}%
!{(\varepsilon_3-\varepsilon_4)(\varepsilon'_3-\varepsilon'_2)-%
!(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_4)})
! \end{align}
! \end{subequations} 
! 
!From the intersections listed in the last equation, the subroutine 
!automatically eliminates those that either do not belong to the surface
!of the tetrahedron or are indefinite.
! 
!We use a parameter 'ndabc' here, since in the later section, for some
!cases we need to get the region of $e1<ef$ and $e2>ef$ by minusing the 
!$e1<ef$ and $e2<ef$ region from the $e1<ef$ region. So, for $ndabc=1$ and
!$ndabc=3$, we get all the nodes for the intersection of the six planes. 
!For $ndabc=2$, we just get that of the five planes except $e2=ef$, to get 
!the region of $e1<ef$

! !USES:

      use polyhedron
      
      implicit none

! !LOCAL VARIABLES:
 
      integer(4) :: i,j,k,inod,iint
      
      real(8) :: edif, fdif, denom, denom1, denom2, denom3, efdif, nt(3),sumnt,efef
      
      real(8), parameter :: zerotol=1.0d-6
      
! !REVISION HISTORY:
!      
! Created 1st. April 2004 by RGA
! Last Revision: 14th. Dec 2004 by XZL

!EOP
!BOC      
!
!     inline functions
!
      edif(i,j)=e(i)-e(j)      
      fdif(i,j)=f(i)-f(j)
      efdif(i,j,k)=edif(i,j)*fdif(k,j)-edif(k,j)*fdif(i,j)
      efef(i,j)=(ef-e(i))*(f(j)-f(i))-(ef-f(i))*(e(j)-e(i))
      intnodes(1:3,1:20) = 0.0d0
      ntype(1)=7
      ntype(2)=11
      ntype(3)=13
      ntype(4)=14
      ntype(5)=19
      ntype(6)=21
      ntype(7)=22
      ntype(8)=25
      ntype(9)=26
      ntype(10)=28
      ntype(11)=35
      ntype(12)=37
      ntype(13)=38
      ntype(14)=41
      ntype(15)=42
      ntype(16)=44
      ntype(17)=49
      ntype(18)=50
      ntype(19)=52
      ntype(20)=56
!
!     The four corners of the tetrahedron:
!     intnodes(1) = plane1^plane2^plane3
!     intnodes(2) = plane1^plane2^plane4
!     intnodes(3) = plane1^plane3^plane4
!     intnodes(4) = plane2^plane3^plane4
!      
      inod=0
      iint=0
      if(ndabc.ne.4) then     
       do i=2,4
         intnodes(5-i,i) = 1.0d0
       enddo
       inod=4 
       iint=4
      else
       continue
      endif

!
!     Intersections of the plane approximating the Fermi surf. at k 
!     with the stright lines defining the edges of the
!     tetrahedron

!      
!     intnodes(5) = plane1^plane2^plane5
!     intnodes(6) = plane1^plane3^plane5
!     intnodes(7) = plane3^plane3^plane5
! 

       do j=2,4
        i=6-j
        iint=iint+1
        denom=edif(i,1)
        nt(1:3)=0.0d0
        if(dabs(denom).gt.zerotol)then
          nt(i-1) = (ef-e(1))/denom
          if((nt(i-1).gt.0.0d0).and.(nt(i-1).le.1.0d0))then
            inod=inod+1
            intnodes(1:3,inod)=nt(1:3)
            ntype(inod)=ntype(iint)
          endif
        endif
       enddo
!      
!     intnodes(8) = plane1^plane4^plane5
!     intnodes(9) = plane2^plane4^plane5
!     intnodes(10) = plane3^plane4^plane5
! 
        do k=1,3
         i=mod(k,3)+1
         j=mod(i,3)+1
         iint=iint+1
         denom=edif(i+1,j+1)
         nt(1:3)=0.0d0
         if(dabs(denom).gt.zerotol)then
          nt(i) = (ef-e(j+1))/denom
          nt(j) = (e(i+1)-ef)/denom
          if((nt(i).ge.0.0d0).and.(nt(j).ge.0.0d0))then
            inod=inod+1
            intnodes(1:3,inod)=nt(1:3)
            ntype(inod)=ntype(iint)
          endif
        endif
       enddo
     
      if((ndabc.ne.2).and.(ndabc.ne.4)) then
!
!     Intersections of the plane approximating the Fermi surf. at k - q
!     with the stright lines defining the edges of the
!     tetrahedron
!
!      write(17,*)
!      
!     intnodes(11) = plane1^plane2^plane6
!     intnodes(12) = plane1^plane3^plane6
!     intnodes(13) = plane3^plane3^plane6
! 
       do j=2,4
        i=6-j
        iint=iint+1
        denom=fdif(i,1)
        nt(1:3)=0.0d0
        if(dabs(denom).gt.zerotol)then
          nt(i-1) = (ef-f(1))/denom
          if((nt(i-1).ge.0.0d0).and.(nt(i-1).le.1.0d0))then
            inod=inod+1
            intnodes(1:3,inod)=nt(1:3)
            ntype(inod)=ntype(iint)
          endif
        endif
       enddo
          
!     intnodes(14) = plane1^plane4^plane6
!     intnodes(15) = plane2^plane4^plane6
!     intnodes(16) = plane3^plane4^plane6
! 
   

       do k=1,3
        i=mod(k,3)+1
        iint=iint+1
        j=mod(i,3)+1
        denom=fdif(i+1,j+1)
        nt(1:3)=0.0d0
        if(dabs(denom).gt.zerotol)then
          nt(i) = (ef-f(j+1))/denom
          nt(j) = (f(i+1)-ef)/denom
          if((nt(i).ge.0.0d0).and.(nt(j).ge.0.0d0))then
            inod=inod+1
            intnodes(1:3,inod)=nt(1:3)
            ntype(inod)=ntype(iint)
          endif
        endif
       enddo
      else
       continue
      endif

      if((ndabc.ne.2).and.(ndabc.ne.4)) then
 
!
!     Intersections of the two planes approximating the Fermi surf. and the 
!     planes defining the tetrahedron surface
!
!      
!     intnodes(17) = plane1^plane5^plane6
!     intnodes(18) = plane2^plane5^plane6
!     intnodes(19) = plane3^plane5^plane6
! 
       do k=1,3
        i=mod(k,3)+1
        iint=iint+1
        j=mod(i,3)+1
        denom1=efdif(i+1,1,j+1)
        denom2=efdif(j+1,1,i+1)
        nt(1:3)=0.0d0
        if((dabs(denom1).gt.zerotol).and.(dabs(denom2).gt.zerotol))then
          nt(i)=efef(1,j+1)/denom1
          nt(j)=efef(1,i+1)/denom2
          sumnt=nt(i)+nt(j)
          if((nt(i).gt.0.0d0).and.(nt(j).gt.0.0d0)                    &
     &       .and.(sumnt.le.1.0d0))then
            inod=inod+1
            intnodes(1:3,inod)=nt(1:3)
            ntype(inod)=ntype(iint)
          endif
        endif
       enddo
!
!     intnodes(20) = plane4^plane5^plane6
! 
      denom1=efdif(2,4,3)
      denom2=efdif(3,4,2)
      denom3=efdif(4,2,3)
      iint=iint+1
       if(((dabs(denom1).gt.zerotol).and.(dabs(denom2).gt.zerotol)).and.(dabs(denom3).gt.zerotol))then
         nt(1)=efef(4,3)/denom1
         nt(2)=efef(4,2)/denom2
         nt(3)=efef(2,3)/denom3
         if((nt(1).gt.0.0d0).and.(nt(2).gt.0.0d0)                      &
     &       .and.(nt(3).gt.0.0d0))then
          inod=inod+1
          intnodes(1:3,inod)=nt(1:3)
          ntype(inod)=ntype(iint)
        endif
       endif
      else
       continue
      endif
      nnod=inod  
 
      end subroutine surfnodes
      
!EOC
