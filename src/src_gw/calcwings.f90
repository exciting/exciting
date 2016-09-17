
subroutine calcwings(ik,iq,iomstart,iomend,ndim,mstart,mend)

    use modmain,   only : evalsv, evalcr, pi, zzero, zone, idxas
    use modgw
    
    implicit none
    ! input variables
    integer, intent(in) :: ik
    integer, intent(in) :: iq
    integer, intent(in) :: iomstart, iomend
    integer, intent(in) :: ndim
    integer, intent(in) :: mstart, mend
    ! local variables
    integer :: nmdim
    integer :: ie1, ie2, ie12
    integer :: is, ia, ic, icg, ias
    integer :: iop, iom, ikp, jk
    real(8) :: edif
    complex(8) :: coefw
    complex(8), allocatable :: pm(:), tvec(:), tmat(:,:)
    
    ! wings prefactor
    coefw = cmplx(-dsqrt(4.d0*pi*vi)*occmax,0.d0,8)
    
    ! irreducible k-point index
    ikp = kset%ik2ikp(ik)
    ! k-q point
    jk = kqset%kqid(ik,iq)
    
    nmdim = ndim*(mend-mstart+1)
    allocate(pm(nmdim))
    allocate(tmat(1:mbsiz,1:nmdim))
    allocate(tvec(1:mbsiz))
    
    !-------------------------------
    ! loop over vector components
    !-------------------------------
    do iop = 1, 3

      !==============================
      ! Sum over states
      !==============================
      ie12 = 0
      do ie2 = mstart, mend
        do ie1 = 1, ndim
          ie12 = ie12+1
          if (ie1<=nomax) then
            !==============================
            ! valence-valence contribution
            !==============================
            edif = evalsv(ie1,ikp)-evalsv(ie2,ikp)
            if (dabs(edif)>1.d-6) then
              pm(ie12) = pmatvv(ie1,ie2,iop)/edif
            else
              pm(ie12) = zzero
            end if
          else
            !==============================
            ! core-valence contribution
            !==============================              
            icg = ie1-nomax
            is = corind(icg,1)
            ia = corind(icg,2)
            ic = corind(icg,3)
            ias = idxas(ia,is)
            edif = evalcr(ic,ias)-evalsv(ie2,ikp)
            if (dabs(edif)>1.d-6) then
              pm(ie12) = pmatcv(icg,ie2,iop)/edif
            else
              pm(ie12) = zzero
            end if
          end if ! core/valence
        end do ! ie1
      end do ! ie2
      
      !==================
      ! Frequency loop
      !==================
      do iom = iomstart, iomend
        !---------
        ! Wing 1
        !---------
        ie12 = 0
        do ie2 = mstart, mend
          do ie1 = 1, ndim
            ie12 = ie12+1
            tmat(1:mbsiz,ie12) = fnm(ie1,ie2,iom,ik)* &
            &                    minmmat(1:mbsiz,ie1,ie2)
          end do ! ie1
        end do ! ie2
        call zgemv('n',mbsiz,nmdim,coefw,tmat,mbsiz,pm,1,zzero,tvec,1)
        epsw1(:,iom,iop) = epsw1(:,iom,iop)+tvec(:)
        !---------
        ! Wing 2
        !---------
        if (freq%fconv=='refreq') then
          ie12 = 0
          do ie2 = mstart, mend
            do ie1 = 1, ndim
              ie12 = ie12+1
              tmat(1:mbsiz,ie12) = conjg(fnm(ie1,ie2,iom,ik))* &
              &                    minmmat(1:mbsiz,ie1,ie2)
            end do ! ie1
          end do ! ie2
          call zgemv('n',mbsiz,nmdim,coefw,tmat,mbsiz,pm,1,zzero,tvec,1)
        end if
        epsw2(:,iom,iop) = epsw2(:,iom,iop)+conjg(tvec(:))
      end do ! iom
      
    end do ! iop
    
    deallocate(tmat)
    deallocate(tvec)
    deallocate(pm)

    return
end subroutine
