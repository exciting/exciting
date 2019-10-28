
subroutine get_evec_gw(vkl, vgkl, evec)
    use modinput
    use modmain, only : ngkmax, nmatmax, &
                        nstfv, nstsv, nspinor, zzero, ldapu
    implicit none
    real(8), intent(in)     :: vkl(3)
    real(8), intent(in)     :: vgkl(3,ngkmax)
    complex(8), intent(out) :: evec(nmatmax,nstfv)

    ! local variables
    integer :: ist, ispn, i, j
    complex(8), allocatable :: evecfvt(:,:)
    complex(8), allocatable :: evecsvt(:,:)

    evec(:,:) = zzero

    if (ldapu /= 0) then
        allocate(evecfvt(nmatmax,nstfv))
        evecfvt = zzero
        allocate(evecsvt(nstsv,nstsv))
        evecsvt = zzero
        call getevecfv(vkl, vgkl, evecfvt)
        call getevecsv(vkl, evecsvt)
        do j = 1, nstsv
            i = 0
            do ispn = 1, nspinor
                do ist = 1, nstfv
                    i = i+1
                    call zaxpy(nmatmax, evecsvt(i,j), &
                               evecfvt(:,ist), 1, &
                               evec(:,j), 1)
                end do
            end do
        end do
        deallocate(evecfvt)
        deallocate(evecsvt)
    else
        call getevecfv(vkl, vgkl, evec)
    end if

    return
end subroutine
