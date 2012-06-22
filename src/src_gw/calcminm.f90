!BOP
!
! !ROUTINE: calcminm
!
! !INTERFACE:
      subroutine calcminm(ik,iq,flag)

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$ 
!of equation \ref{Mdef}.
!
!
! !USES:

      use modinput
      use modmain
      use modgw


! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ik   ! the index of the first k-point
      integer(4), intent(in) :: iq   ! the index of the q-point
      integer(4), intent(in) :: flag ! >0 - transform M^i_{nm} to another basis
      
! !LOCAL VARIABLES:
      integer(4) :: jk      ! the index of the second k-point
      integer(4) :: ikp,jkp
      integer(4) :: bl      ! Angular momentum quantum number of the 
!                             mixed atomic functions (big l)
      integer(4) :: bm      ! z-component angular momentum quantum number
!                             of the mixed atomic function (big m)
      integer(4) :: i
      integer(4) :: ia 
      integer(4) :: ias 
      integer(4) :: igk1,igk2 ! Counter: run over the (L)APW basis 
!                             functions (excluding core states).
      integer(4) :: ilo1 
      integer(4) :: ilo2 
      integer(4) :: imix    ! Counter: runs over all MT-sphere mixed basis 
!                             functions (of all the atoms)
      integer(4) :: io1 
      integer(4) :: io2 
      integer(4) :: irm     ! Counter: runs over the radial mixed basis 
!                             functions of each atom.
      integer(4) :: is 
      integer(4) :: ist1,ist2 ! Counter: run over the eigenvalues in 
!                             the corresponding k-point.
      integer(4) :: l1      ! l-quantum number of the (L)APW eigenfunction 
!                             at k
      integer(4) :: l2      ! l-quantum number of the (L)APW eigenfunction 
!                             at k' = k - q
      integer(4) :: m1      ! m-quantum number of the (L)APW eigenfunction 
!                             at k
      integer(4) :: m2      ! m-quantum number of the (L)APW eigenfunction 
!                             at k' = k - q
      integer(4) :: l2min   ! Lowest allowed value of lambda = | bl - l2 |
      integer(4) :: l2max   ! Highest allowed value of lambda = bl + l2
      integer(4) :: l1m1    ! Counter: Runs over MTS (L)APW basis 
!                             functions (l1,m1)
      integer(4) :: l2m2    ! Counter: Runs over MTS (L)APW basis 
!                               functions l2,m2)
      integer(4) :: igq    ! Counter: Runs over IPW's

      integer(4), allocatable :: igqk12(:,:)   ! Counter: Runs over IPW's

      integer(4), dimension(3):: ikv,ig0        ! Indexes of G_1+G'-G
      
      integer(4) :: ngk1,ngk2
      
      integer(4) :: nmdim
      
      real(8) :: tstart,tend,sqvi,x,arg
      real(8) :: qvec(3)

      real(8) :: angint     ! Angular integral calculated by fourylm
 
      complex(8) :: phs
      complex(8) :: sumterms
      complex(8) :: naterm(apwordmax)
      complex(8) :: nloterm(nlomax)
      
      complex(8), allocatable :: minm(:,:,:)
      complex(8), allocatable :: zzk(:,:)
      complex(8), allocatable :: tmat(:,:),tmat2(:,:),mnn(:,:)
!
! !EXTERNAL ROUTINES: 
      
      real(8), external :: gaunt
      real(8), external :: getcgcoef
      external          :: zgemm

! !INTRINSIC ROUTINES: 

      intrinsic abs
      intrinsic cpu_time
      intrinsic min
      intrinsic nint

     
! !REVISION HISTORY:
! 
! Created  23th. Feb. 2004 by RGA
! Last modified 20th. Jul. 2004 by RGA
! Revisited 29.04.2011 by DIN
!
! !TO DO:
!
! - Check complex and real cases
!
!EOP
!
!BOC
      call cpu_time(tstart)

      allocate(minm(matsiz,nstfv,nstfv))
      minm=zzero
      
      jk=kqid(ik,iq)
      do i=1,3
        qvec(i)=vql(i,iq)
        x=vklnr(i,ik)-vklnr(i,jk)-qvec(i)
        ig0(i)=nint(x)
      enddo
      !write(*,*)'calcminm: ik, jk, ig0: ', ik, jk, ig0
      
      ikp=indkp(ik)
      jkp=indkp(jk)
!
!     -----------------------------------------------
!     Valence and conduction states for both n and n'
!     -----------------------------------------------
!
      imix = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias=idxas(ia,is)
          arg=atposl(1,ia,is)*qvec(1)+atposl(2,ia,is)*qvec(2)+          &
     &        atposl(3,ia,is)*qvec(3)
          phs=cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)
!
!         Loop over mixed functions:
!
          do irm = 1, nmix(ias)
            bl=bigl(ias,irm)

            do bm=-bl,bl
              imix = imix + 1

              do l1=0,input%groundstate%lmaxapw
                l2min=abs(bl-l1)
                l2max=min(bl+l1,input%groundstate%lmaxapw)
                
                do l2=l2min,l2max
                  
                  do m1=-l1,l1
                    l1m1 = l1*l1+l1+m1+1
                    m2=-bm+m1
                    
                    if(abs(m2).le.l2) then
                      l2m2=l2*l2+l2+m2+1
!
!                     Calculate the angular integral:
!
                      !angint=getcgcoef(l2,bl,l1,m2,bm) ! original
                      angint=gaunt(l1,l2,bl,m1,m2,bm)
                      
                      if(abs(angint).gt.1.0d-8) then
!
!                       Loop over eigenfunctions at k
!
                        do ist1 = 1, nstfv

                          naterm=zzero
                          do io2=1,apword(l2,is)
                            do io1=1,apword(l1,is)
                              naterm(io2)= naterm(io2)+                 &
     &                          bradketa(ias,irm,l1,io1,l2,io2,2)*      &
     &                          eveckalm(ist1,io1,l1m1,ias,1)
                            enddo ! io2
                            do ilo1=1,nlorb(is)
                              if(lorbl(ilo1,is).eq.l1)then
                                igk1=ngk(1,ikp)+idxlo(l1m1,ilo1,ias)
                                naterm(io2)=naterm(io2)+                &
     &                            bradketlo(ias,irm,ilo1,l2,io2,2)*     &
     &                            eveck(igk1,ist1,1)
                               endif
                            enddo ! ilo1   
                          enddo ! io2 
                          
                          nloterm=zzero 
                          do ilo2=1,nlorb(is)
                            if(lorbl(ilo2,is).eq.l2)then
                              do io1=1,apword(l1,is)
                                nloterm(ilo2)= nloterm(ilo2)+           &
     &                            bradketa(ias,irm,l1,io1,ilo2,1,3)*    &
     &                            eveckalm(ist1,io1,l1m1,ias,1)
                              enddo ! io1
                              do ilo1=1,nlorb(is)
                                if(lorbl(ilo1,is).eq.l1)then
                                  igk1=ngk(1,ikp)+idxlo(l1m1,ilo1,ias)
                                  nloterm(ilo2)=nloterm(ilo2)+            &
     &                              bradketlo(ias,irm,ilo1,ilo2,1,3)*     &
     &                              eveck(igk1,ist1,1)
                                endif
                              enddo ! ilo1   
                            endif  
                          enddo ! io2  
!
!                         Loop over basis functions at k-q
! 
                          do ist2 = 1, nstfv
                            sumterms=zzero
                            do io2=1,apword(l2,is)
                              sumterms=sumterms+naterm(io2)*eveckpalm(ist2,io2,l2m2,ias,1)
                            enddo ! io2 
                            do ilo2=1,nlorb(is)
                              if(lorbl(ilo2,is).eq.l2)then
                                igk2=ngk(1,jkp)+idxlo(l2m2,ilo2,ias)
                                sumterms=sumterms+nloterm(ilo2)*eveckp(igk2,ist2,1)
                              endif
                            enddo ! ilo2 
                            minm(imix,ist1,ist2) = minm(imix,ist1,ist2)+ &
                           &                       phs*angint*sumterms
                          enddo ! ist2
                        
                        enddo !ist1
                      endif
                    endif
                  enddo ! m1     
                enddo ! l2
              enddo ! l1
            enddo ! bm
          enddo ! irm
        enddo ! ia
      enddo ! is
!
!     --------------------
!     Interstitial region:
!     --------------------
!
!     Loop over the mixed basis functions:
!     
      sqvi=sqrt(vi)
      ngk1=ngknr(1,ik)
      ngk2=ngknr(1,jk)
      allocate(igqk12(ngk1,ngk2))
      allocate(tmat(ngk2,ngk1))
      allocate(tmat2(ngk2,nstfv))
      allocate(mnn(nstfv,nstfv))
      allocate(zzk(1:ngkmax,1:nstfv))

      igqk12(:,:)=0
      do igk1 = 1, ngk1 ! loop over G
        do igk2 = 1, ngk2 ! loop over G'
          ikv(1:3) = ivg(1:3,igkignr(igk1,1,ik)) - &
         &           ivg(1:3,igkignr(igk2,1,jk)) + ig0(1:3)
          if((ikv(1).ge.intgv(1,1)).and.(ikv(1).le.intgv(1,2)).and.       &
         &   (ikv(2).ge.intgv(2,1)).and.(ikv(2).le.intgv(2,2)).and.       &
         &   (ikv(3).ge.intgv(3,1)).and.(ikv(3).le.intgv(3,2)))        then
              igqk12(igk1,igk2)=igigqb(ivgig(ikv(1),ikv(2),ikv(3)),iq)
          endif    
        enddo ! igk2
      enddo ! igk1

      do igq = 1, ngq(iq)
        
        do igk1 = 1, ngk1 ! loop over G
          do igk2 = 1, ngk2 ! loop over G'
            if(igqk12(igk1,igk2).gt.0) then
              tmat(igk2,igk1)=mpwipw(igq,igqk12(igk1,igk2))
            else 
              tmat(igk2,igk1)=zzero
            endif    
          enddo ! igk2
        enddo ! igk1    
        
        zzk(1:ngk1,1:nstfv)=eveck(1:ngk1,1:nstfv,1)
        call zgemm('n','n',ngk2,nstfv,ngk1,zone,tmat,         &
     &             ngk2,zzk,ngkmax,zzero,tmat2,ngk2)

        zzk(1:ngk2,1:nstfv)=eveckp(1:ngk2,1:nstfv,1)
        call zgemm('t','n',nstfv,nstfv,ngk2,zone,zzk,          &
     &             ngkmax,tmat2,ngk2,zzero,mnn,nstfv)

        do ist2 = 1, nstfv  
          do ist1 = 1, nstfv
            minm(igq+locmatsiz,ist1,ist2)=sqvi*mnn(ist2,ist1)
          enddo ! ist2
        enddo ! ist1
      enddo ! igq

      deallocate(igqk12)
      deallocate(tmat)
      deallocate(tmat2)
      deallocate(mnn)
      deallocate(zzk)

!------------------------------------------------------------------    
!     Transform M^i_{nm} to the eigenvectors of the coulomb matrix
!------------------------------------------------------------------

      if(allocated(minmmat))deallocate(minmmat)
      
      nmdim=nstfv*nstfv
      if(flag>0)then
        allocate(minmmat(mbsiz,nstfv,nstfv))
        call zgemm('c','n',mbsiz,nmdim,matsiz, &
       &  zone,barcvm,matsiz,minm,matsiz,zzero,minmmat,mbsiz)
      else
        allocate(minmmat(matsiz,nstfv,nstfv))
        minmmat=minm
      end if ! flag
      
      deallocate(minm)
     
      call cpu_time(tend)
      if (tend.lt.0.0d0)   write(fgw,*) 'warning, tend < 0'
      if (tstart.lt.0.0d0) write(fgw,*) 'warning, tstart < 0'

      return
      end subroutine calcminm
!EOC      


