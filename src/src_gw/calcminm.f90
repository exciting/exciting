!BOP
!
! !ROUTINE: calcminm
!
! !INTERFACE:
      subroutine calcminm(ik,iq)

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$. 
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
      complex(8) :: apwterm(apwordmax)
      complex(8) :: loterm(nlomax)
      
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
      
      if (allocated(minmmat)) deallocate(minmmat)
      allocate(minmmat(matsiz,nstfv,nstfv))
      minmmat = zzero
      
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
          arg=atposl(1,ia,is)*qvec(1)+ &
         &    atposl(2,ia,is)*qvec(2)+ &
         &    atposl(3,ia,is)*qvec(3)
          phs=cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)
!
!         Loop over mixed functions:
!
          do irm = 1, nmix(ias)

            bl=bigl(ias,irm)
            do bm = -bl, bl
              
              imix = imix+1

              do l1=0,input%groundstate%lmaxapw
                
                l2min=abs(bl-l1)
                l2max=min(bl+l1,input%groundstate%lmaxapw)
                do l2=l2min,l2max
                  
                  do m1=-l1,l1
                    l1m1=idxlm(l1,m1)
                    m2=-bm+m1
                    if (abs(m2).le.l2) then
                      l2m2=idxlm(l2,m2)
!
!                     Calculate the angular integral:
!
                      angint=gaunt(l1,l2,bl,m1,m2,bm)
                      
                      if (abs(angint).gt.1.0d-8) then
!
!                       Loop over eigenfunctions at k-q
!
                        do ist2 = 1, nstfv

                          apwterm=zzero
                          do io1=1,apword(l1,is)
                            
                            do io2=1,apword(l2,is)
                              apwterm(io1)= apwterm(io1)+                &
                             &  eveckpalm(ist2,io2,l2m2,ias,1)*          &
                             &  bradketa(ias,irm,l1,io1,l2,io2,2)
                            end do ! io2
                            
                            do ilo2=1,nlorb(is)
                              if (lorbl(ilo2,is).eq.l2) then
                                igk2=ngk(1,jkp)+idxlo(l2m2,ilo2,ias)
                                apwterm(io1)=apwterm(io1)+               &
                               &  eveckp(igk2,ist2,1)*                   &
                               &  bradketa(ias,irm,l1,io1,ilo2,1,3)
                               end if
                            end do ! ilo2
                          
                          end do ! io1 
                          
                          loterm=zzero
                          do ilo1=1,nlorb(is)
                            if (lorbl(ilo1,is).eq.l1) then

                              do io2=1,apword(l2,is)
                                loterm(ilo1)= loterm(ilo1)+           &
                               &  eveckpalm(ist2,io2,l2m2,ias,1)*       &      
                               &  bradketlo(ias,irm,ilo1,l2,io2,2)
                              end do ! io1
                              
                              do ilo2=1,nlorb(is)
                                if (lorbl(ilo2,is).eq.l2) then
                                  igk2=ngk(1,jkp)+idxlo(l2m2,ilo2,ias)
                                  loterm(ilo1)=loterm(ilo1)+            &
                                 &  eveckp(igk2,ist2,1)*                  &
                                 &  bradketlo(ias,irm,ilo1,ilo2,1,3)
                                end if
                              end do ! ilo2   

                            end if
                          end do ! ilo1
!
!                         Loop over basis functions at k
! 
                          do ist1 = 1, nstfv
                            
                            sumterms=zzero
                            
                            do io1=1,apword(l1,is)
                              sumterms=sumterms+eveckalm(ist1,io1,l1m1,ias,1)*apwterm(io1)
                            enddo ! io1
                            
                            do ilo1=1,nlorb(is)
                              if (lorbl(ilo1,is).eq.l1) then
                                igk1=ngk(1,ikp)+idxlo(l1m1,ilo1,ias)
                                sumterms=sumterms+eveck(igk1,ist1,1)*loterm(ilo1)
                              endif
                            enddo ! ilo2 
                            
                            minmmat(imix,ist1,ist2) = minmmat(imix,ist1,ist2) &
                           &                        + phs*angint*sumterms
                          
                          enddo ! ist1
                        
                        enddo !ist2
                      
                      endif ! angint
                    
                    endif ! m2
                    
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
            minmmat(igq+locmatsiz,ist1,ist2)=sqvi*mnn(ist2,ist1)
          enddo ! ist2
        enddo ! ist1
      enddo ! igq

      deallocate(igqk12)
      deallocate(tmat)
      deallocate(tmat2)
      deallocate(mnn)
      deallocate(zzk)
     
      call cpu_time(tend)
      if (tend.lt.0.0d0)   write(fgw,*) 'warning, tend < 0'
      if (tstart.lt.0.0d0) write(fgw,*) 'warning, tstart < 0'

      return
      end subroutine calcminm
!EOC      


