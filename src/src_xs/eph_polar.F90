
subroutine eph_polar(nw, w, head)
    use modmain, only: natmtot, natoms, nspecies, idxas, omega, spmass
    use constants, only: zi
    use constants, only: fourpi
    use m_getunit
    implicit none
    integer(4), intent(in)    :: nw
    complex(8), intent(in)    :: w(nw)
    complex(8), intent(inout) :: head(3,3,nw)
    ! local
    integer(4) :: fid, nmode, imode, i, j, is, ia, ias, js, ja, iw
    character(len=40) :: fname
    real(8), allocatable :: phfreq(:)
    real(8), allocatable :: zborn(:,:,:)
    complex(8), allocatable :: phevec(:,:,:)
    complex(8), allocatable :: zvec(:,:), v(:)
    complex(8) :: elat(3,3,nw), om
    real(8) :: vnorm

    !--------------------------------------
    ! read Born effective charge tensors
    !--------------------------------------
    if (allocated(zborn)) deallocate(zborn)
    allocate(zborn(3,3,natmtot))
    zborn(:,:,:) = 0.d0

    fname = 'BORN.OUT'
    call getunit(fid)
    open(unit=fid, file=trim(fname), status='old', action='read')
    read(fid,*) ! skip comment
    do is = 1, nspecies
        do ia = 1, natoms(is)
            ias = idxas(ia,is)
            read(fid,*) js, ja, zborn(:,:,ias)
            ! write(*,*) is, ia, zborn(:,:,ias)
        end do
    end do

    !------------------------------------------
    ! read phonon frequencies and eigenvectors
    !------------------------------------------
    nmode = 3*natmtot
    if (allocated(phfreq)) deallocate(phfreq)
    allocate(phfreq(nmode))
    phfreq(:) = 0.d0

    if (allocated(phevec)) deallocate(phevec)
    allocate(phevec(3,natmtot,nmode))
    phevec(:,:,:) = 0.d0

    fname = 'PHONON.OUT'
    call getunit(fid)
    open(unit=fid, file=trim(fname), status='old', action='read')
    read(fid,*) ! q-point info
    read(fid,*) ! skip empty line

    do imode = 1, nmode
        read(fid,*) i, phfreq(imode)
        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia,is)
                do i = 1, 3
                    ! read eigenvectors
                    read(fid,*) js, ja, j, phevec(i,ias,imode)
                    phevec(i,ias,imode) = phevec(i,ias,imode) / sqrt(spmass(is))
                end do
            end do
        end do
        read(fid,*) ! skip empty line
    end do

    !-------------------------------------------------------------------------------
    ! Compute the e-ph contribution to the head of the dielectric function
    ! according to Eqs.(54-55) of X. Gonze and C. Lee, PRB 55, 10355 (1997)
    !-------------------------------------------------------------------------------
    if (allocated(zvec)) deallocate(zvec)
    allocate(zvec(3,nmode))
    zvec(:,:) = 0.d0

    do imode = 1, nmode
        do i = 1, 3
            ! sum over atoms
            do ias = 1, natmtot
                zvec(i,imode) = zvec(i,imode) + dot_product(zborn(i,:,ias), phevec(:,ias,imode))
            end do
        end do
        ! write(*,'(a,i4,3f12.4)') 'zvec=', imode, dble(zvec(:,imode))
    end do

    elat(:,:,:) = 0.d0
    do iw = 1, nw
        do i = 1, 3
            do j = 1, 3
                ! sum over phonon modes
                do imode = 1, nmode
                    if (abs(phfreq(imode)) > 1.d-4) then ! skip acoustic modes
                        ! elat(i,j,iw) = elat(i,j,iw) + conjg(zvec(i,imode))*zvec(j,imode) / (phfreq(imode)**2 - w(iw)**2)
                        elat(i,j,iw) = elat(i,j,iw) + conjg(zvec(i,imode))*zvec(j,imode) / phfreq(imode)**2
                    end if
                end do
            end do
        end do
    end do
    elat = elat * fourpi/omega

    deallocate(zvec)
    deallocate(zborn)
    deallocate(phfreq)
    deallocate(phevec)

    ! write(*,*)
    ! do i = 1, 3
    !     write(*,'(2f16.4,4x,2f16.4,4x,2f16.4)') (elat(i,j,1), j = 1, 3)
    ! end do
    ! write(*,*)
    ! stop 'eph_polar'

    head = head + elat

end subroutine
