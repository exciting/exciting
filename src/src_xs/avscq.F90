! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: acscq
! !INTERFACE:
subroutine avscq(iqr, n, nmax, scrn, scieff)
! !USES:
  use modinput, only: input
  use mod_qpoint, only: iqmap
  use modxs, only: ivqr
  use modmpi
  use invert
! !INPUT/OUTPUT PARAMETERS:
!   IN:
!   iqr, integer             : Index of the reduced q-point to be considered
!   n, integer               : Number of G+q vectors for current q-point
!   nmax, integer            : Maximum number of G+q vectors over all q-points
!   scrn(n,n), complex(8)    : Body of dielectric matrix
!   OUT:
!   scieff(nmax,nmax), complex(8) : Screened coulomb potential 
!   
! !DESCRIPTION:
!   By default this just calculates
!   $\frac{\tilde{\varepsilon}_{G,G'}(q,\omega=0)}{|G+q||G'+q|}$.
!   There are options controlling how the dielectric matrix should be 
!   inverted. There are also possible averaging and extrapolation schemes. 
!
! !REVISION HISTORY:
! Added to documentation scheme. (Aurich)
!
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: iqr, n, nmax
  complex(8), intent(in) :: scrn(n, n)
  complex(8), intent(out) :: scieff(nmax, nmax)

  ! Local variables
  integer :: iqrnr, j1, j2, flg_analytic, flg
  real(8) :: clwt
  logical :: doavrage

  ! Find reduced q-point in non-reduced set
  iqrnr = iqmap(ivqr(1, iqr), ivqr(2, iqr), ivqr(3, iqr))

  ! Invert dielectric matrix
  call zinvert_hermitian(input%xs%bse%scrherm, scrn, scieff(:n, :n))

  select case(trim(input%xs%bse%cuttype))
    case("none")
      flg_analytic = 0
    case("0d")
      flg_analytic = 4
    case("2d")
      flg_analytic = 5
    case default
      write(*,*) "Error(avscq): Invalid cuttype"
      call terminate
  end select

  doavrage = input%xs%bse%sciavqhd .or.&
           & input%xs%bse%sciavqwg .or.&
           & input%xs%bse%sciavqbd

  if(flg_analytic /= 0  .and. doavrage) then
    write(*,*) "Error(avscq): Cutoff currently only&
      & supported for analytic expression."
    call terminate
  end if

  do j1 = 1, n
    do j2 = 1, j1

        ! Average head of Coulomb potential for q != 0 
        ! (Default: sciavqhd = false)
      if((input%xs%bse%sciavqhd .and. (j1 .eq. 1) .and. (j2 .eq. 1))&
        ! Average wings of Coulomb potential for q != 0 
        ! (Default: sciavqwg = false)
        & .or. (input%xs%bse%sciavqwg .and. (j1 .ne. 1) .and. (j2 .eq. 1))&
        & .or. (input%xs%bse%sciavqwg .and. (j1 .eq. 1) .and. (j2 .ne. 1))&
        ! Average body of Coulomb potential for q != 0 
        ! (Default: sciavqbd = false)
        & .or. (input%xs%bse%sciavqbd .and. (j1 .ne. 1) .and. (j2 .ne. 1))) then
        ! Numerical averaging on grids with extrapolation to continuum
        flg = 2
      else
        ! Default case:
        ! Analytic expression, no averaging
        flg = flg_analytic
      end if

      ! Generate the (averaged) symmetrized coulomb potential
      ! flag = 0 : clwt = v^{1/2}(G,q) * v^{1/2}(G',q)
      call genwiqggp(flg, iqrnr, j1, j2, clwt)

      ! Multiply with (averaged) coulomb potential
      scieff(j1, j2) = scieff(j1, j2) * clwt

      ! Make explicitly hermitian.
      scieff(j2, j1) = conjg(scieff(j1, j2))

    end do
  end do

end subroutine avscq
!EOC
