!
    subroutine calcmicm(ik,iq)
!
    use modinput
    use modmain
    use modgw

    implicit none
    integer(4), intent(in)  :: ik   ! k-point
    integer(4), intent(in)  :: iq   ! q-point

    integer(4) :: jk
    integer(4) :: bl, bm
    integer(4) :: ia, is, ias
    integer(4) :: igk2
    integer(4) :: io2, ilo2
    integer(4) :: icg, ic, ie2
    integer(4) :: imix, irm, im
    integer(4) :: l1, m1, l2, m2, l2m2 
    integer(4) :: l2min, l2max
     
    real(8) :: arg, kvec(3)
    real(8) :: tstart,tend

    complex(8) :: phs, angint, sum

    external :: zgemm
    real(8), external :: gaunt

    call timesec(tstart)

    !---------------------------------------------------
    ! n => core state, n' => Valence or conduction state
    !---------------------------------------------------
    
    if (allocated(micmmat)) deallocate(micmmat)
    allocate(micmmat(locmatsiz,ncg,nstfv))
    micmmat = zzero
    
    jk = kqid(ik,iq)
    kvec(:) = vklnr(:,jk)

    ! loop over core states
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(icg,is,ia,ic,l1,m1,ias,arg,phs,im,irm,bl,bm,imix,l2min,l2max,l2,m2,l2m2,angint,ie2,io2,sum,ilo2,igk2)
!$OMP DO
#endif 
    do icg = 1, ncg
      is  = corind(icg,1)
      ia  = corind(icg,2)
      ic  = corind(icg,3)
      l1  = corind(icg,4)
      m1  = corind(icg,5)      
      
      ias = idxas(ia,is)
      arg = atposl(1,ia,is)*kvec(1)+ &
      &     atposl(2,ia,is)*kvec(2)+ &
      &     atposl(3,ia,is)*kvec(3)
      phs = cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)
          
      ! Loop over mixed functions:
      im = 0
      do irm = 1, nmix(ias)
        bl = bigl(ias,irm)
        do bm = -bl, bl
          im = im+1
          imix = locmixind(ias,im)

          ! sum over l2m2
          l2min = iabs(bl-l1)
          l2max = min(bl+l1,input%groundstate%lmaxapw)
          do l2 = l2min, l2max
            m2 = -bm+m1
            !m2 = bm+m1 !<-- FHI-gap
            if (iabs(m2) <= l2) then
              l2m2 = idxlm(l2,m2)
                                  
              ! Angular integral
              angint = gaunt(l1,bl,l2,m1,bm,m2)
              !angint = getgauntcoef(bl,l2,l1,bm,m2)
              !angint = getgauntcoef(l1,bl,l2,m1,bm) !<-- FHI-gap
                
              if (abs(angint)>1.0d-8) then
                
                ! loop over valence states
                do ie2 = 1, nstfv
                  sum = zzero
                  !------------------
                  ! APW contribution
                  !------------------
                  do io2 = 1, apword(l2,is)
                    sum = sum + &
                    &     eveckpalm(ie2,io2,l2m2,ias,1)* &
                    &     bradketc(ias,irm,ic,l2,io2,2)
                  end do
                  !------------------
                  ! LO contribution
                  !------------------
                  do ilo2 = 1, nlorb(is)
                    if (lorbl(ilo2,is)==l2) then
                      igk2 = ngknr(1,jk)+idxlo(l2m2,ilo2,ias)
                      sum = sum + &
                      &     eveckp(igk2,ie2,1)* &
                      &     bradketc(ias,irm,ic,ilo2,1,3)
                    end if
                  end do  
                  
                  micmmat(imix,icg,ie2) = micmmat(imix,icg,ie2)+phs*angint*sum
                
                end do ! ie2

              end if ! angint
              
            end if  
          end do ! l2
          
        end do ! bm  
      end do ! irm
        
    end do ! icg
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    

    !------------------------------------------------------------------
    ! timing
    call timesec(tend)

    return
end subroutine
