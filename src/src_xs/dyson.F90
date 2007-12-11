
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dyson
  implicit none
contains
  
  subroutine dyson(n,s0,k,s)
    !
    ! solve Dyson's equation S = S0 + S0*K*S
    ! for S by inversion.
    !
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: s0(:,:), k(:,:)
    complex(8), intent(out) :: s(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'dyson'
    complex(8),parameter :: zone=(1.d0,0.d0),zzero=(0.d0,0.d0)
    complex(8), allocatable :: mt(:,:),zwork(:)
    integer, allocatable :: ipiv(:)
    integer :: shs0(2),shk(2),shs(2),nmin,nmax,j,lwork,info

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
       stop
    end if

    ! allocate
    allocate(mt(n,n))
    allocate(ipiv(n))
    lwork=2*n
    allocate(zwork(lwork))

    ! calculate matrix -S0*K
    call zgemm('n','n', n, n, n, -zone, s0, n, k, n, zzero, mt, n )

    ! calculate matrix T=[1 - S0*K]
    forall(j=1:n) mt(j,j)=mt(j,j)+1.d0
    
    ! invert matrix T
    call zgetrf(n,n,mt,n,ipiv,info)
    if (info.ne.0) then
       write(*,*)
       write(*,'("Error(",a,"): zgetrf returned non-zero info : ",I8)') &
            thisnam,info
       write(*,*)
       stop
    end if
    call zgetri(n,mt,n,ipiv,zwork,lwork,info)
    if (info.ne.0) then
       write(*,*)
       write(*,'("Error(",a,"): zgetri returned non-zero info : ",I8)') &
            thisnam,info
       write(*,*)
       stop
    end if
    
    ! calculate matrix S=T^-1*S0
    call zgemm('n','n', n, n, n, zone, mt, n, s0, n, zzero, s, n )

    deallocate(mt,ipiv,zwork)

  end subroutine dyson
  
end module m_dyson
