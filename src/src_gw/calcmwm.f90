
subroutine calcmwm(nstart, nend, mstart, mend, minm)

    use modinput
    use constants, only: zone, zzero, pi
    use modgw,   only: vi, kqset, Gamma, singc1, singc2, mbsiz, &
    &                  epsilon, epsh, epsw1, epsw2, freq, mwm
    implicit none

    ! input variables
    integer(4), intent(in) :: nstart, nend
    integer(4), intent(in) :: mstart, mend
    complex(8), intent(in) :: minm(mbsiz, nstart:nend, mstart:mend)

    ! local variables
    integer(4) :: iom
    integer(4) :: ie1, ie2
    real(8)    :: wkq
    real(8)    :: vi4pi, coefs1, coefs2, coefs3
    complex(8), allocatable :: wm(:)
    complex(8), external    :: zdotu, zdotc
    external zhemm

    wkq    = 1.d0/dble(kqset%nkpt)
    
    if (Gamma) then
      vi4pi  = 4.d0*pi*vi
      coefs1 = singc1*sqrt(vi4pi)
      coefs2 = singc2*vi4pi
    end if 

    !-------------------------------------------------
    ! calculate \sum_{ij} M^i_{nm}* W^c_{ij} M^j_{nm}
    !-------------------------------------------------
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(wm,iom,ie1,ie2)
#endif
    allocate(wm(mbsiz))
    do iom = 1, freq%nomeg
      do ie2 = mstart, mend
#ifdef USEOMP
!$OMP DO SCHEDULE(DYNAMIC)
#endif
        do ie1 = nstart, nend
          call zhemv( 'u', mbsiz, zone, epsilon(:,:,iom), mbsiz, &
          &           minm(:,ie1,ie2), 1, zzero, wm, 1)
          mwm(ie1,ie2,iom) = wkq*zdotc(mbsiz,minm(:,ie1,ie2),1,wm,1)
          if ((Gamma).and.(ie1==ie2)) then
            mwm(ie1,ie2,iom) = mwm(ie1,ie2,iom) + &
            &                  coefs2*epsh(1,1,iom) + &
            &                  coefs1*(zdotu(mbsiz, minm(:,ie1,ie2), 1, epsw2(:,1,iom), 1) + &
            &                          zdotc(mbsiz, minm(:,ie1,ie2), 1, epsw1(:,1,iom), 1))
          end if ! singular term
        end do
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif
      end do
    end do ! iom
    deallocate(wm)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

    return
end subroutine
