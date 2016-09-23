!BOP
!
!!ROUTINE: calcmicm
!
!!INTERFACE:
!
    subroutine calcmicm(ik,iq,cstart,cend,mstart,mend,micm)
!
! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{cm}(\vec{k},\vec{q})$ 
!(c stands for core-state)
!
!!USES:
    use modinput
    use modmain,           only : pi, idxas, idxlm, idxlo, nlorb, apword, &
    &                             lorbl, zzero
    use modgw,             only : kqset, Gkset, fdebug, time_micm
    use mod_bands,         only : eveckp, eveckpalm
    use mod_product_basis, only : locmatsiz, nmix, bigl, locmixind, bradketc
    use mod_core_states,   only : corind, ncg
    use mod_misc_gw,       only : atposl
    use mod_gaunt_coefficients

!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in)  :: ik   ! k-point
    integer(4), intent(in)  :: iq   ! q-point
    integer(4), intent(in)  :: cstart, cend  ! range of core states
    integer(4), intent(in)  :: mstart, mend  ! range of m states
    complex(8), intent(out) :: micm(locmatsiz,ncg,mstart:mend)

!!LOCAL VARIABLES:
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

!!EXTERNAL ROUTINES:
    external :: zgemm
    real(8), external :: gaunt

!!REVISION HISTORY:
! 
! Created  23th. Feb. 2004 by RGA
! Last modified 20th. July 2004 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC

    call timesec(tstart)

    !---------------------------------------------------
    ! n => core state, n' => Valence or conduction state
    !---------------------------------------------------
    
    micm(:,:,:) = zzero
    
    jk = kqset%kqid(ik,iq)
    kvec(:) = kqset%vkl(:,jk)

    ! loop over core states
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(icg,is,ia,ic,l1,m1,ias,arg,phs,im,irm,bl,bm,imix,l2min,l2max,l2,m2,l2m2,angint,ie2,io2,sum,ilo2,igk2)
!$OMP DO
#endif 
    do icg = cstart, cend
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
        bl = bigl(irm,ias)
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
              !angint = gaunt(l1,bl,l2,m1,bm,m2)
              angint = getgauntcoef(bl,l2,l1,bm,m2)
              !angint = getgauntcoef(l1,bl,l2,m1,bm) !<-- FHI-gap
                
              if (abs(angint)>1.0d-8) then
                
                ! loop over valence states
                do ie2 = mstart, mend
                  sum = zzero
                  !------------------
                  ! APW contribution
                  !------------------
                  do io2 = 1, apword(l2,is)
                    sum = sum + &
                    &     eveckpalm(ie2,io2,l2m2,ias)* &
                    &     bradketc(2,irm,ic,l2,io2,ias)
                  end do
                  !------------------
                  ! LO contribution
                  !------------------
                  do ilo2 = 1, nlorb(is)
                    if (lorbl(ilo2,is)==l2) then
                      igk2 = Gkset%ngk(1,jk)+idxlo(l2m2,ilo2,ias)
                      sum = sum + &
                      &     eveckp(igk2,ie2)* &
                      &     bradketc(3,irm,ic,ilo2,1,ias)
                    end if
                  end do  
                  
                  micm(imix,icg,ie2) = micm(imix,icg,ie2)+phs*angint*sum
                
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
    time_micm = time_micm+tend-tstart

    return
end subroutine
!EOC      


