

! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dysonsym
  implicit none
contains
  
!BOP
! !ROUTINE: dysonsym
! !INTERFACE:
  subroutine dysonsym(n,s0,k,s)
! !USES:
    use invert
! !INPUT/OUTPUT PARAMETERS:
!   n     : matrix size of local field effects (in,integer)
!   s0    : S0 matrix (in,complex(:,:))
!   k     : kernel matrix multiplied by S0 from both sides (in,complex(:,:))
!   s     : S (solution) matrix (in,complex(:,:))
! !DESCRIPTION:
!   Solve symmetric form of Dyson's equation
!     $$   S = S_0 + S_0 (1 + S0^-1 K S0^-1) S  $$
!   for $S$ by inversion;
!     $$ S = S_0\left[ S_0(1-S0) - T\right]^{-1} S_0. $$
!   The inversion is carried out using the LAPACK routines {\tt zgetrf} and
!   {\tt zgetri}.
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: s0(:,:), k(:,:)
    complex(8), intent(out) :: s(:,:)
    ! local variables
    character(*), parameter :: thisnam='dysonsym'
    complex(8),parameter :: zone=(1.d0,0.d0),zzero=(0.d0,0.d0)
    complex(8), allocatable :: mt(:,:),mt2(:,:)
    integer :: shs0(2),shk(2),shs(2),nmin,nmax,j

    ! check matrix sizes
    shs0=shape(s0)
    shk=shape(k)
    shs=shape(s)
    nmin=minval((/shs0,shk,shs/))
    nmax=maxval((/shs0,shk,shs/))
    if ((nmin.ne.nmax).or.(nmin.lt.n)) then
       write(*,'("Error(",a,"): inconsistent matrix sizes")') trim(thisnam)
       write(*,'("  n :",i9)') n
       write(*,'("  S0:",2i9)') shs0
       write(*,'("  K :",2i9)') shk
       write(*,'("  S :",2i9)') shs
       call terminate
    end if

    ! allocate
    allocate(mt(n,n),mt2(n,n))

    ! calculate matrix 1-S0
    mt(:,:)=zzero
    forall(j=1:n) mt(j,j)=1.d0
    mt(:,:)=mt(:,:)-s0(:,:)
    
    ! calculate S0(1-S0)
    call zgemm('n','n', n, n, n, zone, s0, n, mt, n, zzero, mt2, n )
   
    ! calculate S0(1-S0) - K
    mt2(:,:)=mt2(:,:)-k(:,:)
    
    ! calculate [S0(1-S0) - K]^-1 =: Y
    call zinvert_lapack(mt2,mt)

    ! calculate S0 Y
    call zgemm('n','n', n, n, n, zone, s0, n, mt, n, zzero, mt2, n )
    
    ! calculate solution S = S0 Y S0 = S0 [S0(1-S0) - K]^(-1) S0
    call zgemm('n','n', n, n, n, zone, mt2, n, s0, n, zzero, s, n )
    
    deallocate(mt,mt2)

  end subroutine dysonsym
!EOC

end module m_dysonsym



