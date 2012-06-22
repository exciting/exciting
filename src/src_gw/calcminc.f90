!BOP
!
! !ROUTINE: calcmicm
!
! !INTERFACE:
      subroutine calcminc(ik,iq,flag)

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{cm}(\vec{k},\vec{q})$ 
!of equation \ref{Mdef}. c stands for core-state
!
!
! !USES:
      
      use modinput
      use modmain
      use modgw

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ik  ! k-point
      integer(4), intent(in) :: iq  ! q-point
      integer(4), intent(in) :: flag ! >0 - transform M^i_{cm} to another basis

! !LOCAL VARIABLES:

      integer(4) :: jk    ! the k-point index
      integer(4) :: ikp

      integer(4) :: bl    ! Angular momentum quantum number of the 
!                           atomic mixed function irm.(big l)
      integer(4) :: bm    ! z-component angular momentum quantum number
!                           of the atomic function irm. (big m)
      integer(4) :: i
      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: ias0
      integer(4) :: igk2
      integer(4) :: ilo2
      integer(4) :: io2
      integer(4) :: is
      integer(4) :: icore ! Counter: runs over the core states of a 
!                           given atom.
      integer(4) :: icg   ! Runs over all core states of all atoms
      integer(4) :: ist   ! Counter: run over the eigenvalues in 
!                           the corresponding k-point.
      integer(4) :: imix  ! Counter: runs over all MT-sphere mixed basis 
!                           functions (of all the atoms)
      integer(4) :: irm   ! Counter: runs over the radial mixed basis 
!                           functions of each atom.
      integer(4) :: l1    ! l-quantum number of the (L)APW eigenfunction 
!                           at k
      integer(4) :: l2    ! l-quantum number of the (L)APW eigenfunction 
!                           at k' = k - q
      integer(4) :: l1m1,l2m2   ! Counter: Runs over MTS (L)APW basis 
!                                 functions (l1,m1 and l2,m2 
!                                 respectively)
      integer(4) :: l2min,l2max ! minumum and maximum (respectively)
!                                 possible values of l2       
      integer(4) :: m1   ! m-quantum number of the (L)APW eigenfunction 
!                          at k
      integer(4) :: m2   ! m-quantum number of the (L)APW eigenfunction 
!                          at k' = k - q
      integer(4) :: im
      
      integer(4) :: ncdim
      
      real(8) :: arg
      real(8) :: tstart,tend
      real(8), dimension(3) :: kqvec
      
      complex(8) :: phs
      complex(8) :: angint  ! Angular integral calculated by fourylm
      complex(8) :: suml2m2 ! intermediate sum (over (L)APW basis 
!                             functions [l2, m2])
      complex(8) :: sumterms
      complex(8) :: apwterm  
      complex(8) :: loterm  
      
      complex(8), allocatable :: minc(:,:,:)
!
! !EXTERNAL ROUTINES: 
      
      real(8), external :: gaunt
      real(8), external :: getcgcoef

! !INTRINSIC ROUTINES: 

      intrinsic iabs
      intrinsic isign
      intrinsic min
!
! !REVISION HISTORY:
! 
! Created  23th. Feb. 2004 by RGA
! Last modified 20th. July 2004 by RGA
! Revisited: May 2011 by DIN
!
! !TO DO:
!
! - Check complex and real cases
!
!EOP
!
!BOC
!
!     ---------------------------------------------------
!      n => core state, n' => Valence or conduction state
!     ---------------------------------------------------

      call cpu_time(tstart)
      
      allocate(minc(locmatsiz,1:nstfv,ncg))
      minc=zzero

      jk=kqid(ik,iq)
      do i=1,3
        kqvec(i)=-vklnr(i,ik)
      enddo
      ikp=indkp(ik)
!
!     Loop over core states
!
      ias0=0
      do icg = 1, ncg
        is=corind(icg,1)
        ia=corind(icg,2)
        ias=idxas(ia,is)
    
        arg=atposl(1,ia,is)*kqvec(1)+atposl(2,ia,is)*kqvec(2)+        &
     &      atposl(3,ia,is)*kqvec(3)
        phs=cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)

!       reset istlms if it is a new atom 
        if(ias.ne.ias0)then
          l1m1=0
          ias0=ias
        endif  
        icore=corind(icg,3)
        l1=corind(icg,4)
        m1=corind(icg,5)
        l1m1=l1m1+1
!
!       Loop over APW states
!
        do ist = 1, nstfv
!
!         Loop over mixed functions:
! 
          imix=0
          do irm = 1, nmix(ias)
            bl=bigl(ias,irm)
            do bm=-bl,bl
              imix=imix+1
              im=locmixind(ias,imix)
              
              l2min=iabs(bl-l1)
              l2max=min(bl+l1,input%groundstate%lmaxapw)
              suml2m2 = zzero

              do l2=l2min,l2max
                m2=bm+m1
                if(iabs(m2).le.l2)then
                  l2m2=l2*l2+l2+m2+1
                  
                  ! Gaunt coefficient
                  !angint=getcgcoef(l2,bl,l1,m2,bm) ! original
                  angint=gaunt(l2,l1,bl,m2,m1,bm)
                  
                  if(abs(angint).gt.1.0d-8)then
                    apwterm=zzero
                    do io2=1,apword(l2,is)
                      apwterm = apwterm + &
     &                          bradketc(ias,irm,icore,l2,io2,2)*  &
     &                          eveckalm(ist,io2,l2m2,ias,1)
                    enddo
                    loterm=zzero
                    do ilo2=1,nlorb(is)
                      if(lorbl(ilo2,is).eq.l2)then
                        igk2=ngk(1,ikp)+idxlo(l2m2,ilo2,ias)
                        loterm = loterm + &
     &                           bradketc(ias,irm,icore,ilo2,1,3)* &
     &                           eveck(igk2,ist,1)
                      endif
                    enddo  
                    sumterms = apwterm + loterm
                    suml2m2 = suml2m2 + angint * sumterms
                  endif
                endif  
              enddo ! l2

              minc(im,ist,icg) = suml2m2*phs

            enddo ! bm  
          enddo ! irm

        enddo ! ist

      enddo ! icg
      
!------------------------------------------------------------------    
!     Transform M^i_{nc} to the eigenvectors of the coulomb matrix
!------------------------------------------------------------------
      
      if(allocated(mincmat))deallocate(mincmat)

      ncdim=nstfv*ncg
      if(flag>0)then
        allocate(mincmat(mbsiz,1:nstfv,ncg))
        call zgemm('c','n',mbsiz,ncdim,locmatsiz, &
       &  zone,barcvm,matsiz,minc,locmatsiz,zzero,mincmat,mbsiz)
      else
        allocate(mincmat(locmatsiz,1:nstfv,ncg))
        mincmat=minc
      end if ! flag
      
      deallocate(minc)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'

      end subroutine calcminc
!EOC      


