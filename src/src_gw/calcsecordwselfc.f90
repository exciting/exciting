subroutine calcsecordwselfc(ikp,iq,mdim)

    use modgw
    use mod_mpi_gw
    implicit none
    
    integer, intent(in) :: ikp
    integer, intent(in) :: iq
    integer, intent(in) :: mdim
    
    integer :: iom
    integer :: ie1, ie2, nmdim
    integer :: ik, jk, jkp
    integer :: ia, is, ias, ic, icg
    real(8) :: vi4pi, coefs1, coefs2
    real(8) :: wkq, enk
    complex(8) :: sc
    complex(8) :: xnm(1:freq%nomeg)
    complex(8) :: wm(mbsiz,ibgw:nbgw,1:mdim)
    
    complex(8), external :: freqconv
    complex(8), external :: zdotu, zdotc
    external zhemm
    
    vi4pi = 4.d0*pi*vi
    coefs1 = singc1*sqrt(vi4pi)
    coefs2 = singc2*vi4pi
    wkq = 1.d0/dble(kqset%nkpt)
    
    do iom = 1, freq%nomeg
      ! calculate \sum_{j} W^c_{ij}M^j_{nm}
      nmdim = (nbgw-ibgw+1)*mdim
      call zhemm('l','u',mbsiz,nmdim, &
      &          zone,vPv(:,:,iom),mbsiz,minmmat,mbsiz, &
      &          zzero,wm,mbsiz)
      do ie2 = 1, mdim
        do ie1 = ibgw, nbgw
          mwm(ie1,ie2,iom) = wkq*zdotc(mbsiz,minmmat(:,ie1,ie2),1,wm(:,ie1,ie2),1)
          if ((Gamma).and.(ie1==ie2)) then
            mwm(ie1,ie2,iom) = mwm(ie1,ie2,iom) + &
            &  coefs2*vPvh(iom) + &
            &  coefs1*(zdotu(mbsiz,minmmat(:,ie1,ie2),1,vPvw2(:,iom),1) + &
            &          zdotc(mbsiz,minmmat(:,ie1,ie2),1,vPvw1(:,iom),1))
          end if ! singular term
        end do
      end do
    end do ! iom
    
    !===========================
    ! Frequency convolution
    !===========================
    ! k-point
    ik = kset%ikp2ik(ikp)
    ! k-q point
    jk = kqset%kqid(ik,iq)
    jkp = kset%ik2ikp(jk)

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iom,ie1,sc,ie2,enk,xnm,icg,is,ia,ic,ias)
!$OMP DO
#endif    
    do iom = 1, freq%nomeg
      ! loop over states
      do ie1 = ibgw, nbgw
      
        ! sum over states
        sc = zzero
        do ie2 = 1, mdim
        
          if (ie2<=nstsv) then
            !============================= 
            ! Valence electron contribution
            !============================= 
            enk = evalsv(ie2,jkp)-efermi
            xnm(1:freq%nomeg) = mwm(ie1,ie2,1:freq%nomeg)
          else
            !============================= 
            ! Core electron contribution
            !=============================
            icg = ie2-nstsv
            is = corind(icg,1)
            ia = corind(icg,2)
            ic = corind(icg,3)
            ias = idxas(ia,is)
            enk = evalcr(ic,ias)-efermi
            xnm(1:freq%nomeg) = mwm(ie1,ie2,1:freq%nomeg)
          end if ! val/cor
          
          !------------------------
          ! Frequency convolution
          !------------------------
          sc = sc+freqconv(iom,freq%nomeg,freq%freqs(iom), &
          &                enk,xnm,freq%freqs,freq%womeg)
                  
        end do ! ie2
    
        selfecw2(ie1,iom,ikp) = selfecw2(ie1,iom,ikp)+sc
    
      end do ! ie1
    end do ! iom
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif        
    
    return
end subroutine
