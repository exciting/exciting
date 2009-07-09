


! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dyson
  implicit none
contains

!BOP
! !ROUTINE: dyson
! !INTERFACE:


subroutine dyson(n, s0, k, s)
! !USES:
    use invert
! !INPUT/OUTPUT PARAMETERS:
!   n     : matrix size of local field effects (in,integer)
!   s0    : S0 matrix (in,complex(:,:))
!   k     : kernel matrix (in,complex(:,:))
!   s     : S (solution) matrix (in,complex(:,:))
! !DESCRIPTION:
!   Solve Dyson's equation
!     $$   S = S_0 + S_0 K S  $$
!   for $S$ by inversion;
!     $$ S = \left[ 1 + S_0 K \right]^{-1} S_0. $$
!   The inversion is carried out using the LAPACK routines {\tt zgetrf} and
!   {\tt zgetri}.
!
! !REVISION HISTORY:
!   Created March 2005 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: s0(:, :), k(:, :)
    complex(8), intent(out) :: s(:, :)
    ! local variables
    character(*), parameter :: thisnam='dyson'
    complex(8), parameter :: zone=(1.d0, 0.d0), zzero=(0.d0, 0.d0)
    complex(8), allocatable :: mt(:, :)
    integer :: shs0(2), shk(2), shs(2), nmin, nmax, j

    ! check matrix sizes
    shs0=shape(s0)
    shk=shape(k)
    shs=shape(s)
    nmin=minval((/shs0, shk, shs/))
    nmax=maxval((/shs0, shk, shs/))
    if ((nmin.ne.nmax).or.(nmin.lt.n)) then
       write(*, '("Error(", a, "): inconsistent matrix sizes")') trim(thisnam)
       write(*, '("  n :", i9)') n
       write(*, '("  S0:", 2i9)') shs0
       write(*, '("  K :", 2i9)') shk
       write(*, '("  S :", 2i9)') shs
       call terminate
    end if

    ! allocate
    allocate(mt(n, n))

    ! calculate matrix -S0*K
    call zgemm('n', 'n', n, n, n, -zone, s0, n, k, n, zzero, mt, n )

    ! calculate matrix T=[1 - S0*K]
    forall(j=1:n) mt(j, j)=mt(j, j)+1.d0

    ! invert matrix T
    call zinvert_lapack(mt, mt)

    ! calculate matrix S=T^-1*S0
    call zgemm('n', 'n', n, n, n, zone, mt, n, s0, n, zzero, s, n )

    deallocate(mt)

  end subroutine dyson
!EOC

end module m_dyson
