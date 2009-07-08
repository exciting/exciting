



subroutine orthogonalise(n, nstfv, evecfv, ldv, s)
implicit none
integer, intent(in)::n, nstfv, ldv
complex(8), intent(inout)::evecfv(ldv, nstfv)
complex(8), intent(in)::s(n, nstfv)

integer :: ist, jst, i
complex(8)::zt1
complex(8), external:: zdotc
real(8)::t1
! perform Gram-Schmidt orthonormalisation
  do ist=1, nstfv
    do jst=1, ist-1
      zt1=-zdotc(n, evecfv(1, jst), 1, s(1, ist), 1)
      call zaxpy(n, zt1, evecfv(1, jst), 1, evecfv(1, ist), 1)
    end do
  end do
end subroutine
