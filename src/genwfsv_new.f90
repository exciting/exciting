
subroutine genwfsv_new(ik, ist1, ist2, apwalm, evecfv, evecsv, wfmt, wfir)

    use modinput
    use modmain
    implicit none

    ! arguments
    integer,    intent(in)  :: ik
    integer,    intent(in)  :: ist1, ist2
    complex(8), intent(in)  :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)  :: evecfv(nmatmax,nstfv)
    complex(8), intent(in)  :: evecsv(nstsv,nstsv)
    complex(8), intent(out) :: wfmt(lmmaxapw,nrmtmax,natmtot,nspinor,nstsv)
    complex(8), intent(out) :: wfir(ngrtot,nspinor,nstsv)
    ! local variables
    integer    :: ispn, is, ia, ias
    integer    :: i, j, n, ist, igp
    real(8)    :: t1
    complex(8) :: zt1
    ! allocatable arrays
    logical,    allocatable :: done (:)
    complex(8), allocatable :: wfmt1(:,:)

    wfmt(:,:,:,:,:) = 0.d0
    wfir(:,:,:) = 0.d0

    allocate(done(nstfv))
    allocate (wfmt1(lmmaxapw,nrcmtmax))

!--------------------------------!
!     muffin-tin wavefunction    !
!--------------------------------!
    
    do is = 1, nspecies
      n = lmmaxapw*nrmt(is)
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        done(:) = .false.
        do j = ist1, ist2
          if (input%groundstate%tevecsv) then
            ! generate spinor wavefunction from second-variational eigenvectors
            wfmt(:,:,ias,:,j) = 0.d0
            i = 0
            do ispn = 1, nspinor
              do ist = 1, nstfv
                i = i + 1
                zt1 = evecsv(i,j)
                if (abs(zt1)>1.d-8) then
                  if (.not.done(ist)) then
                    call wavefmt(1, input%groundstate%lmaxapw, is, ia, &
                    &            ngk(1,ik), apwalm, evecfv(:,ist), lmmaxapw, wfmt1)
                    done(ist) = .true.
                  end if
                  ! add to spinor wavefunction
                  call zaxpy(n, zt1, wfmt1, 1, &
                  &          wfmt(:,:,ias,ispn,j), 1)
                end if
              end do ! ist
            end do ! ispn
          else
            ! spin-unpolarised wavefunction
            call wavefmt(1, input%groundstate%lmaxapw, is, ia, &
            &            ngk(1,ik), apwalm, evecfv(:,j), lmmaxapw, &
            &            wfmt(:,:,ias,1,j))
          end if
        end do ! j
      end do ! ia
    end do ! is
    deallocate(done, wfmt1)

!-----------------------------------!
!     interstitial wavefunction     !
!-----------------------------------!
    ! normalization factor is included
    t1 = 1.d0/dsqrt(omega)
    do j = ist1, ist2
      if (input%groundstate%tevecsv) then
        ! generate spinor wavefunction from second-variational eigenvectors
        i = 0
        do ispn = 1, nspinor
          do ist = 1, nstfv
            i = i + 1
            zt1 = evecsv(i,j)
            if (abs(zt1)>1.d-8) then
              zt1 = t1*zt1
              do igp = 1, ngk(1,ik)
                wfir(igp,ispn,j) = wfir(igp,ispn,j)+zt1*evecfv(igp,ist)
              end do
            end if
          end do
        end do
      else
        ! spin-unpolarised wavefunction
        do igp = 1, ngk(1,ik)
          wfir(igp,1,j) = t1*evecfv(igp,j)
        end do
      end if
    end do ! j
    
    return
end subroutine
