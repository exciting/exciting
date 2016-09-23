!BOP
!
!!ROUTINE: calcminm
!
!!INTERFACE:
!
subroutine calcminm(ik,iq,nstart,nend,mstart,mend,minm)
!
! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$ 
!
!!USES:
    use modinput
    use modmain,               only : nspecies, natoms, idxas, idxlm, idxlo, &
    &                                 zzero, zone, intgv, apword, nlorb, lorbl, &
    &                                 nlomax, pi, apwordmax, nmatmax
    use modgw,                 only : kqset, Gkset, Gqbarc, Gqset, Gset, fdebug, time_minm
    use mod_bands,             only : eveck, eveckp, eveckalm, eveckpalm
    use mod_product_basis,     only : nmix, bigl, bradketa, bradketlo, mpwipw, &
    &                                 matsiz, locmatsiz, mbindex
    use mod_misc_gw,           only : vi, atposl
    use mod_gaunt_coefficients

!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: ik   ! the index of the first k-point
    integer(4), intent(in) :: iq   ! the index of the q-point
    integer(4), intent(in) :: nstart, nend  ! range of n states
    integer(4), intent(in) :: mstart, mend  ! range of m states
    complex(8), intent(out):: minm(matsiz,nstart:nend,mstart:mend)
      
!!LOCAL VARIABLES:
    integer(4) :: jk
    integer(4) :: bl, bm
    integer(4) :: i, ia, is, ias 
    integer(4) :: igk1, igk2 
    integer(4) :: io1, io2, ilo1, ilo2 
    integer(4) :: imix, irm
    integer(4) :: ie1, ie2
    integer(4) :: l1, m1, l2, m2, l1m1, l2m2
    integer(4) :: l2min, l2max
    integer(4) :: ig, igq
    integer(4), allocatable :: igqk12(:,:)
    integer(4), dimension(3):: ikv, ig0 ! Indexes of G_1+G'-G
    integer(4) :: ngk1, ngk2
    integer(4) :: ndim, mdim, nmdim
      
    real(8) :: sqvi, x, arg, angint
    real(8) :: qvec(3)
    real(8) :: tstart, tend, tmt, tir
    
    real(8) :: t1, t2
    
    complex(8) :: phs
    complex(8) :: sumterms
    complex(8) :: apwterm(apwordmax)
    complex(8) :: loterm(nlomax)
    
    complex(8), allocatable :: tmat(:,:), tmat2(:,:), mnn(:,:)

!!EXTERNAL ROUTINES: 
    external :: zgemm
    real(8), external :: gaunt

!!REVISION HISTORY:
! 
! Created  23th. Feb. 2004 by RGA
! Last modified 20th. Jul. 2004 by RGA
! Revisited 29.04.2011 by DIN
!
!EOP
!BOC

    call timesec(tstart)
    
    !-----------------------------------------------
    ! Valence and conduction states for both n and n'
    !-----------------------------------------------
    minm(:,:,:) = zzero

    ! index ranges
    ndim = nend-nstart+1
    mdim = mend-mstart+1
    nmdim = ndim*mdim
      
    jk = kqset%kqid(ik,iq)
    do i = 1, 3
      qvec(i) = kqset%vql(i,iq)
      x = kqset%vkl(i,ik)-kqset%vkl(i,jk)-qvec(i)
      ig0(i) = nint(x)
    end do
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(imix,is,ia,ias,irm,bl,bm,arg,phs,l1,m1,l1m1,l2min,l2max,l2,m2,l2m2,angint,ie2,apwterm,io1,io2,ilo2,igk2,loterm,ilo1,ie1,sumterms,igk1)
!$OMP DO
#endif
    ! loop over MT product basis functions
    do imix = 1, locmatsiz
    
      is  = mbindex(imix,1)
      ia  = mbindex(imix,2)
      irm = mbindex(imix,3)
      bl  = mbindex(imix,4)
      bm  = mbindex(imix,5)
      
      ias = idxas(ia,is)
      arg = atposl(1,ia,is)*qvec(1)+ &
      &     atposl(2,ia,is)*qvec(2)+ &
      &     atposl(3,ia,is)*qvec(3)
      phs = cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)

      ! Sum over l1m1 and l2m2
      do l1 = 0, input%groundstate%lmaxapw
        do m1 = -l1, l1
          l1m1 = idxlm(l1,m1)
                
          l2min = abs(bl-l1)
          l2max = min(bl+l1,input%groundstate%lmaxapw)
          do l2 = l2min, l2max
            m2 = -bm+m1
            if (abs(m2).le.l2) then
              l2m2 = idxlm(l2,m2)

              ! Angular integral
              !angint = gaunt(l1,bl,l2,m1,bm,m2)
              angint = getgauntcoef(l2,bl,l1,m2,bm)
              if (abs(angint)<1.d-8) cycle
            
              ! Loop over eigenfunctions at k-q
              do ie2 = mstart, mend
                  
                !------------------
                ! APW contribution
                !------------------
                apwterm(:) = zzero
                do io1 = 1, apword(l1,is)
                      
                  do io2 = 1, apword(l2,is)
                    apwterm(io1) = apwterm(io1)+ &
                    &              eveckpalm(ie2,io2,l2m2,ias)* &
                    &              bradketa(2,irm,l1,io1,l2,io2,ias)
                  end do ! io2
                      
                  do ilo2 = 1, nlorb(is)
                    if (lorbl(ilo2,is)==l2) then
                      igk2 = Gkset%ngk(1,jk)+idxlo(l2m2,ilo2,ias)
                      apwterm(io1) = apwterm(io1)+ &
                      &              eveckp(igk2,ie2)* &
                      &              bradketa(3,irm,l1,io1,ilo2,1,ias)
                    end if
                  end do ! ilo2
                    
                end do ! io1 
                    
                !------------------
                ! LO contribution
                !------------------
                loterm(:) = zzero
                do ilo1 = 1, nlorb(is)
                  if (lorbl(ilo1,is)==l1) then
                  
                    do io2 = 1, apword(l2,is)
                      loterm(ilo1) = loterm(ilo1)+ &
                      &              eveckpalm(ie2,io2,l2m2,ias)* &      
                      &              bradketlo(2,irm,ilo1,l2,io2,ias)
                    end do ! io1
                        
                    do ilo2 = 1, nlorb(is)
                      if (lorbl(ilo2,is)==l2) then
                        igk2 = Gkset%ngk(1,jk)+idxlo(l2m2,ilo2,ias)
                        loterm(ilo1) = loterm(ilo1)+ &
                        &              eveckp(igk2,ie2)* &
                        &              bradketlo(3,irm,ilo1,ilo2,1,ias)
                      end if
                    end do ! ilo2
                  
                  end if
                end do ! ilo1
                  
                ! Loop over basis functions at k
                do ie1 = nstart, nend
                      
                  sumterms = zzero
                  do io1 = 1, apword(l1,is)
                    sumterms = sumterms + &
                    &          eveckalm(ie1,io1,l1m1,ias)*apwterm(io1)
                  end do ! io1
                      
                  do ilo1 = 1, nlorb(is)
                    if (lorbl(ilo1,is)==l1) then
                      igk1 = Gkset%ngk(1,ik)+idxlo(l1m1,ilo1,ias)
                      sumterms = sumterms + &
                      &          eveck(igk1,ie1)*loterm(ilo1)
                    end if
                  end do ! ilo2
                  
                  !--------------------------
                  ! sum up all contributions
                  !--------------------------
                  minm(imix,ie1,ie2) = minm(imix,ie1,ie2)+phs*angint*sumterms
                    
                end do ! ie1
                  
              end do !ie2
                  
            end if ! m2
          end do ! l2
                  
        end do ! m1
      end do ! l1
            
    end do ! imix
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

    call timesec(tmt)
    write(*,*) 'calcminm, mt',tmt-tstart 
    !======================
    ! Interstitial region
    !======================
    
    ! Loop over the mixed basis functions:
    sqvi = sqrt(vi)
    ngk1 = Gkset%ngk(1,ik)
    ngk2 = Gkset%ngk(1,jk)
    
    allocate(igqk12(ngk1,ngk2))
    igqk12(:,:) = 0
    
    do igk1 = 1, ngk1 ! loop over G
      do igk2 = 1, ngk2 ! loop over G'
        ikv(1:3) = Gset%ivg(1:3,Gkset%igkig(igk1,1,ik)) - &
        &          Gset%ivg(1:3,Gkset%igkig(igk2,1,jk)) + ig0(1:3)
        if((ikv(1).ge.Gset%intgv(1,1)).and.(ikv(1).le.Gset%intgv(1,2)).and. &
        &  (ikv(2).ge.Gset%intgv(2,1)).and.(ikv(2).le.Gset%intgv(2,2)).and. &
        &  (ikv(3).ge.Gset%intgv(3,1)).and.(ikv(3).le.Gset%intgv(3,2)))  then
          ig = Gset%ivgig(ikv(1),ikv(2),ikv(3))
          igqk12(igk1,igk2) = Gqbarc%igigk(ig,1,iq)
          !write(*,*) igk1, igk2, ig, igqk12(igk1,igk2)
        end if
      end do ! igk2
    end do ! igk1

!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(igq,igk1,igk2,tmat,tmat2,mnn,ie1,ie2)
!#endif

    allocate(tmat(ngk2,ngk1))
    allocate(tmat2(ngk2,ndim))
    allocate(mnn(mdim,ndim))

!#ifdef USEOMP
!!$OMP DO
!#endif    
    do igq = 1, Gqset%ngk(1,iq)
        
      do igk1 = 1, ngk1 ! loop over G+k
        do igk2 = 1, ngk2 ! loop over G'+k'
          if (igqk12(igk1,igk2) > 0) then
            tmat(igk2,igk1) = mpwipw(igq,igqk12(igk1,igk2))
          else 
            tmat(igk2,igk1) = zzero
          end if    
        end do ! igk2
      end do ! igk1    
      
      call zgemm( 'n', 'n', &
      &           ngk2, ndim, ngk1, &
      &           zone, &
      &           tmat, ngk2, &
      &           eveck(1:ngk1,nstart:nend), ngk1, &
      &           zzero, tmat2, ngk2)

      call zgemm( 't', 'n', &
      &           mdim, ndim, ngk2, &
      &           zone, &
      &           eveckp(1:ngk2,mstart:mend), ngk2, &
      &           tmat2, ngk2, &
      &           zzero, mnn, mdim)

      do ie2 = mstart, mend  
        do ie1 = nstart, nend
          minm(locmatsiz+igq,ie1,ie2) = sqvi*mnn(ie2-mstart+1,ie1-nstart+1)
        end do ! ie2
      end do ! ie1
      
    end do ! igq
!#ifdef USEOMP    
!!$OMP END DO
!#endif

    deallocate(tmat)
    deallocate(tmat2)
    deallocate(mnn)
    
!#ifdef USEOMP
!!$OMP END PARALLEL
!#endif    

    deallocate(igqk12)
    
    !------------------------------------------------------------------
    ! timing
    call timesec(tend)
    time_minm = time_minm+tend-tstart
    write(*,*) 'calcminm, ir', tend-tmt
    !write(*,*) ' minmmat ik, iq: ', ik, iq
    !do imix = 1, matsiz, matsiz/10
    !  do ie2 = 1, 4
    !    do ie1 = 1, 14
    !      write(*,*) imix, ie1, ie2, minm(imix,ie1,ie2)
    !    enddo ! ist2
    !  enddo ! ist1
    !end do
  
    return
end subroutine
!EOC      

