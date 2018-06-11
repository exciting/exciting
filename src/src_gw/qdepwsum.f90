
subroutine qdepwsum(iq,iomstart,iomend,ndim)
      
    use modmain
    use modgw

    implicit none
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: iomstart, iomend
    integer(4), intent(in) :: ndim
      
    integer(4) :: iom, n, m
    integer(4) :: ik, jk, ikp, jkp
    integer(4) :: ia, is, ias, ic, icg
    real(8)    :: de, wkp, e0, eta, ff
    complex(8) :: z1, z2
    real(8)    :: tstart,tend
    complex(8), allocatable :: om(:)

    call timesec(tstart)

    if (allocated(fnm)) deallocate(fnm)
    allocate(fnm(1:ndim,numin:nstsv,iomstart:iomend,1:kqset%nkpt))
    fnm(:,:,:,:) = zzero

    if (allocated(om)) deallocate(om)
    allocate(om(iomstart:iomend))

    select case (freq%fconv)
      case('refreq')
        om(iomstart:iomend) = freq%freqs(iomstart:iomend)
        eta = input%gw%scrcoul%swidth
      case('imfreq')
        om(iomstart:iomend) = zi*freq%freqs(iomstart:iomend)
        eta = 0.d0
      case default
        stop 'Not supported option!'
    end select

    wkp = 1.0d0 / dble(kqset%nkpt)

    do ik = 1, kqset%nkpt
      jk  = kqset%kqid(ik,iq)
      ikp = kset%ik2ikp(ik)
      jkp = kset%ik2ikp(jk)

      do n = 1, ndim

        if (n <= nomax) then
          e0 = evalsv(n,ikp)
        else
          icg = n - nomax
          is  = corind(icg,1)
          ia  = corind(icg,2)
          ic  = corind(icg,3)
          ias = idxas(ia,is)
          e0  = evalcr(ic,ias)
        end if

        do m = numin, nstsv
          
          do iom = iomstart, iomend
            ff = 1.d0 !occsv(n,ikp)/occmax * ( 1.d0 - occsv(m,jkp)/occmax )
            de = evalsv(m,jkp) - e0
            z1 = om(iom) - de + zi*eta
            z2 = om(iom) + de - zi*eta
            fnm(n,m,iom,ik) = ff * (1.d0/z1 - 1.d0/z2) * wkp
          end do ! iom

        end do ! m

      end do ! n
      
    end do ! ik

    deallocate(om)
    
    !-------------------------
    ! Debugging info
    !-------------------------
    if (input%gw%debug) then
      write(fdebug,*)'------------------------------------------------------'
      write(fdebug,*)'       convolution weights for iq =',iq
      write(fdebug,*)'------------------------------------------------------'
      write(fdebug,*)
      iom = 1
      do ik = 1, kqset%nkpt
        do n = 1, ndim
        do m = numin, nstsv
          write(fdebug,10) n, m, iom, ik, fnm(n,m,iom,ik)
        end do
        end do
      end do
      10 format('  n =',i4,' m =',i4,' iom =',i4,' ik =',i4,' fnm =',2g16.8)
    end if ! debug

    call timesec(tend)
    time_bzinit = time_bzinit+tend-tstart

end subroutine
