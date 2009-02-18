

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

  complex(8), allocatable :: solv(:),s0row(:),s0col(:),u(:,:),vh(:,:),work(:)
  real(8), allocatable :: singv(:),rwork(:)
  integer, allocatable :: ipiv(:)
  integer :: info,lwork
  real(8) :: eps
  
  complex(8), external :: zdotu
    
    
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
   
    ! calculate X := S0(1-S0) - K
    mt2(:,:)=mt2(:,:)-k(:,:)
    
!    ! calculate [S0(1-S0) - K]^-1 =: Y = X^-1
!    call zinvert_lapack(mt2,mt)

!    ! calculate S0 Y
!    call zgemm('n','n', n, n, n, zone, s0, n, mt, n, zzero, mt2, n )
    
!    ! calculate solution S = S0 Y S0 = S0 [S0(1-S0) - K]^(-1) S0
!    call zgemm('n','n', n, n, n, zone, mt2, n, s0, n, zzero, s, n )
    
    
    
goto 100    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    ! F. Sottile, PhD thesis, p. 167 (Appendix E)
    ! solve linear system of equations instead of direct inversion
    
    ! first column of S0 is RHS of system of equations
    allocate(solv(n),s0row(n),s0col(n))
    s0row(:)=s0(1,:)
    s0col(:)=s0(:,1)
    solv(:)=s0col(:)
    
    !------------------------------------ solve linear system of equations
    allocate(ipiv(n))
    call zgetrf(n,n,mt2,n,ipiv,info)
    if (info.ne.0) then
       write(*,*)
       write(*,'("Error(",a,"): zgetrf returned non-zero info : ",I8)') &
            thisnam,info
       write(*,*)
       call terminate
    end if
    call zgetrs('n', n, 1, mt2, n, ipiv, solv, n, info )
    if (info.ne.0) then
       write(*,*)
       write(*,'("Error(",a,"): zgetrs returned non-zero info : ",I8)') &
            thisnam,info
       write(*,*)
       call terminate
    end if
    deallocate(ipiv)
    !------------------------------------
    
    ! calculate S_00 = S0 * Solv
    s(1,1)=zdotu(n,solv,1,s0row,1)
    deallocate(solv,s0row,s0col)
100 continue
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    lwork=3*n
    allocate(singv(n),u(n,n),vh(n,n),work(lwork),rwork(5*n))

    ! try SVD for inversion of matrix X = S0(1 - S0) - K
    call ZGESVD( 'a', 'a', n, n, mt2, n, singv, u, n, vh, n, work, lwork, rwork, info )
    if (info.ne.0) then
       write(*,*)
       write(*,'("Error(",a,"): zgesvd returned non-zero info : ",I8)') &
            thisnam,info
       write(*,*)
       call terminate
    end if
   
   ! invert singular values above cutoff
   eps=1.d-8
    do j=1,n
   	if (singv(j).lt.eps) then
	   singv(j)=0.d0
	else
	   singv(j)=1.d0/singv(j)
	end if
	! multiply singular values with U-matrix
	u(:,j)=u(:,j)*singv(j)
    end do
    ! inverse of matrix X^+:
    call zgemm('n','n', n, n, n, zone, u, n, vh, n, zzero, s, n )
    s=conjg(transpose(s))
    
    ! left and right multiply with S0
    ! calculate S0 Y
    call zgemm('n','n', n, n, n, zone, s0, n, s, n, zzero, mt2, n )
    
    ! calculate solution S = S0 Y S0 = S0 [S0(1-S0) - K]^(-1) S0
    call zgemm('n','n', n, n, n, zone, mt2, n, s0, n, zzero, s, n )


    deallocate(mt,mt2,singv,u,vh,work,rwork)

  end subroutine dysonsym
!EOC

end module m_dysonsym



