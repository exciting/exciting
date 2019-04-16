
subroutine get_evec_gw(vkl, vgkl, evec)
    use modinput
    use modmain, only : ngkmax_ptr, nmatmax_ptr, nspnfv, &
                        nstfv, nstsv, nspinor, zzero, ldapu
    implicit none
    real(8), intent(in)     :: vkl(3)
    real(8), intent(in)     :: vgkl(3,ngkmax_ptr,nspnfv)
    complex(8), intent(out) :: evec(nmatmax_ptr,nstfv,nspnfv)

    ! local variables
    integer :: ist, ispn, i, j
    complex(8), allocatable :: evecfvt(:,:)
    complex(8), allocatable :: evecsvt(:,:)

    evec = zzero
    if (ldapu /= 0) then

        allocate(evecfvt(nmatmax_ptr,nstfv))
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
                    call zaxpy(nmatmax_ptr, evecsvt(i,j), &
                               evecfvt(:,ist), 1, &
                               evec(:,j,ispn), 1)
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
