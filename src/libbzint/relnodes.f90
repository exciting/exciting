
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: relnodes
!
! !INTERFACE:
      subroutine relnodes
!     
! !DESCRIPTION:
!
! This subroutine select the relevant intnodes, that is, only those defining
! the polyhedron that has to be integrated. The intnodes selected are those for
! which $\varepsilon(\vec{k})\le \varepsilon_F$ and 
! $\varepsilon'(\vec{k})\ge \varepsilon_F$ if ndabc=1. $\varepsilon(\vec{k})
! \le \varepsilon_F$ when ndabc=2.  $\varepsilon(\vec{k})\le \varepsilon_F$ and
! $\varepsilon'(\vec{k})\le \varepsilon_F$ when ndabc=3. So the first region equals
! the second region minused by the third region.

!
! !USES:

      use polyhedron

      implicit none


! !LOCAL VARIABLES:

      integer(4) :: inod, jnod, maxnod
      
      real(8) :: e1, e2, zerotol
      
      real(8), dimension(3) :: ncor
      
      real(8), external :: tlinap

! !REVISION HISTORY:
!
! Created 22nd. April 2004 by RGA
! Last revised 14th. Dec 2004 by XZL
!
!EOP
!BOC
      zerotol=1.00d-15 
      inod=1
      maxnod=nnod
      select case(ndabc)
      case(1)
!-----------------------------------------------------------------------------------------------
! in this case, we select the nodes to get the region when e1<efer and e2>efer
!-----------------------------------------------------------------------------------------------
       inod=1
       maxnod=nnod
       do while (inod.le.maxnod)
        ncor(1:3)=intnodes(1:3,inod)
        e1=tlinap(ncor,e)
        e2=tlinap(ncor,f)
        if(((e1.gt.ef).and.(dabs(e1-ef).gt.zerotol)).or.((e2.lt.ef).and.(dabs(e2-ef).gt.zerotol)))then
          do jnod=inod+1,maxnod
            intnodes(1:3,jnod-1)=intnodes(1:3,jnod)
            ntype(jnod-1)=ntype(jnod)
          enddo
          maxnod=maxnod-1
        else
          inod=inod+1
        endif
       enddo
       do inod=maxnod+1,nnod
        intnodes(1:3,inod)=0.0d0
        ntype(inod)=0
       enddo

      case(2)
!-----------------------------------------------------------------------------------------------
! in this case, we select the nodes to get the region when e1<efer
!-----------------------------------------------------------------------------------------------
       do while (inod.le.maxnod)
        ncor(1:3)=intnodes(1:3,inod)
        e1=tlinap(ncor,e)
        e2=tlinap(ncor,f)
        if((e1.gt.ef).and.(dabs(e1-ef).gt.zerotol))then
          do jnod=inod+1,maxnod
            intnodes(1:3,jnod-1)=intnodes(1:3,jnod)
            ntype(jnod-1)=ntype(jnod)
          enddo
          maxnod=maxnod-1
        else
          inod=inod+1
        endif
       enddo     
       do inod=maxnod+1,nnod
        intnodes(1:3,inod)=0.0d0
        ntype(inod)=0
       enddo

      case(3)
!-----------------------------------------------------------------------------------------------
! In this case, we get the region when e1<efer and e1<efer
!-----------------------------------------------------------------------------------------------
       inod=1
       maxnod=nnod
       do while (inod.le.maxnod)
        ncor(1:3)=intnodes(1:3,inod)
        e1=tlinap(ncor,e)
        e2=tlinap(ncor,f)
        if(((e1.gt.ef).and.(dabs(e1-ef).gt.zerotol)).or.((e2.gt.ef).and.(dabs(e2-ef).gt.zerotol)))then
          do jnod=inod+1,maxnod
            intnodes(1:3,jnod-1)=intnodes(1:3,jnod)
            ntype(jnod-1)=ntype(jnod)
          enddo
          maxnod=maxnod-1
        else
          inod=inod+1
        endif
       enddo
       do inod=maxnod+1,nnod
        intnodes(1:3,inod)=0.0d0
        ntype(inod)=0
       enddo


      end select
  
      nnod=maxnod

      end subroutine relnodes
      
!EOC 
